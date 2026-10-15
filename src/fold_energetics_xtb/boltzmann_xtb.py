from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence

HARTREE_TO_KCAL_MOL = 627.509474
GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083


@dataclass(frozen=True)
class ConformerState:
    conformer_id: str
    xtb_energy_hartree: float
    solvation_gbsa_kcal_mol: float
    pocket_rmsd_angstrom: float
    hbond_count: int


@dataclass(frozen=True)
class BoltzmannEnsembleResult:
    lowest_energy_conformer: str
    boltzmann_weights: Dict[str, float]
    ensemble_strain_kcal_mol: float
    conformational_entropy_cal_mol_k: float
    free_energy_penalty_kcal_mol: float


class GFN2xTBBoltzmannEvaluator:
    """Computes canonical Boltzmann conformer weights, strain penalties, and PyMOL visualization scripts."""

    def __init__(self, temperature_kelvin: float = 298.15) -> None:
        if temperature_kelvin <= 0.0:
            raise ValueError("temperature_kelvin must be strictly positive")
        self.temperature_kelvin = temperature_kelvin
        self.rt_kcal = GAS_CONSTANT_KCAL_MOL_K * temperature_kelvin

    def evaluate_ensemble(
        self,
        conformers: Sequence[ConformerState],
        bioactive_conformer_id: str | None = None,
    ) -> BoltzmannEnsembleResult:
        if not conformers:
            raise ValueError("conformers cannot be empty")

        total_energies_kcal: List[Tuple[str, float]] = []
        for c in conformers:
            e_kcal = c.xtb_energy_hartree * HARTREE_TO_KCAL_MOL + c.solvation_gbsa_kcal_mol
            total_energies_kcal.append((c.conformer_id, e_kcal))

        min_id, min_e = min(total_energies_kcal, key=lambda item: item[1])
        rel_energies = {cid: e - min_e for cid, e in total_energies_kcal}

        # Log-sum-exp stabilized Boltzmann weights
        neg_beta_e = {cid: -rel_e / self.rt_kcal for cid, rel_e in rel_energies.items()}
        max_log = max(neg_beta_e.values())
        unnorm = {cid: math.exp(val - max_log) for cid, val in neg_beta_e.items()}
        z_part = sum(unnorm.values())
        weights = {cid: round(u / z_part, 6) for cid, u in unnorm.items()}

        # Ensemble-averaged internal strain energy above global minimum
        ensemble_strain = sum((u / z_part) * rel_energies[cid] for cid, u in unnorm.items())

        # Gibbs conformational entropy in cal/(mol*K)
        r_cal = GAS_CONSTANT_KCAL_MOL_K * 1000.0
        entropy = -r_cal * sum(
            (u / z_part) * math.log(max(u / z_part, 1e-15)) for u in unnorm.values()
        )

        target_id = bioactive_conformer_id if bioactive_conformer_id in rel_energies else min_id
        target_prob = max(unnorm[target_id] / z_part, 1e-15)
        free_energy_penalty = -self.rt_kcal * math.log(target_prob)

        return BoltzmannEnsembleResult(
            lowest_energy_conformer=min_id,
            boltzmann_weights=weights,
            ensemble_strain_kcal_mol=round(ensemble_strain, 4),
            conformational_entropy_cal_mol_k=round(entropy, 4),
            free_energy_penalty_kcal_mol=round(free_energy_penalty, 4),
        )

    @staticmethod
    def generate_pymol_pocket_script(
        uniprot_id: str,
        pocket_residues: Sequence[int],
        top_ligand_id: str,
    ) -> str:
        if not pocket_residues:
            raise ValueError("pocket_residues cannot be empty")
        resi_sel = "+".join(str(r) for r in pocket_residues)
        lines = [
            f"# Headless PyMOL Ray-Trace Script for {uniprot_id} ({top_ligand_id})",
            "bg_color white",
            "set ray_trace_mode, 1",
            "set ray_shadows, 0",
            f"fetch AF-{uniprot_id}-F1-model_v4, async=0",
            "hide everything",
            "show cartoon, all",
            "spectrum b, red_yellow_green_cyan_blue, minimum=50, maximum=95",
            f"select catalytic_pocket, resi {resi_sel}",
            "show sticks, catalytic_pocket",
            "color olive, catalytic_pocket",
            "zoom catalytic_pocket, 8",
            f"png {uniprot_id}_{top_ligand_id}_pocket.png, width=1600, height=1200, dpi=300, ray=1",
        ]
        return "\n".join(lines) + "\n"
