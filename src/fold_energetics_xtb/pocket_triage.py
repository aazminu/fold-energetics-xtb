from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import List, Sequence


@dataclass(frozen=True)
class CandidateLigand:
    chembl_id: str
    smiles: str
    pchembl_value: float
    xtb_strain_kcal_mol: float
    pocket_rmsd_angstrom: float
    hbond_count: int


@dataclass(frozen=True)
class PocketEvaluation:
    uniprot_id: str
    pocket_mean_plddt: float
    pae_domain_coherence: float
    druggable_pocket_valid: bool
    top_ligand_id: str
    composite_binding_score: float


class StructuralPocketTriage:
    """Evaluates AlphaFold pocket rigid-domain quality and GFN2-xTB ligand energetics."""

    def __init__(self, min_plddt: float = 70.0, max_pae_angstrom: float = 5.0) -> None:
        self.min_plddt = min_plddt
        self.max_pae_angstrom = max_pae_angstrom

    def score_pocket_and_ligands(
        self,
        uniprot_id: str,
        pocket_plddts: Sequence[float],
        pairwise_pae: Sequence[float],
        ligands: Sequence[CandidateLigand],
    ) -> PocketEvaluation:
        mean_plddt = sum(pocket_plddts) / max(len(pocket_plddts), 1)
        coherent_pairs = sum(1 for p in pairwise_pae if p <= self.max_pae_angstrom)
        coherence = coherent_pairs / max(len(pairwise_pae), 1)
        valid_pocket = mean_plddt >= self.min_plddt and coherence >= 0.75

        best_id = "NONE"
        best_score = -1e9
        for lig in ligands:
            # Penalize conformer strain > 2.5 kcal/mol and high pocket RMSD
            strain_penalty = max(0.0, lig.xtb_strain_kcal_mol - 1.5) * 0.65
            geom_bonus = math.exp(-lig.pocket_rmsd_angstrom / 2.0) + 0.18 * lig.hbond_count
            score = (
                1.25 * lig.pchembl_value
                + 1.40 * geom_bonus
                - strain_penalty
                + (0.50 if valid_pocket else -2.0)
            )
            if score > best_score:
                best_score = score
                best_id = lig.chembl_id

        return PocketEvaluation(
            uniprot_id=uniprot_id,
            pocket_mean_plddt=round(mean_plddt, 2),
            pae_domain_coherence=round(coherence, 4),
            druggable_pocket_valid=valid_pocket,
            top_ligand_id=best_id,
            composite_binding_score=round(best_score, 4),
        )


def summarize_pocket_ensemble(evaluations: Sequence[PocketEvaluation]) -> dict[str, object]:
    """Summarize multi-pocket druggability fraction, mean pLDDT, and top composite binding score."""
    if not evaluations:
        return {
            "total_pockets": 0,
            "druggable_count": 0,
            "druggable_fraction": 0.0,
            "mean_plddt": 0.0,
            "best_uniprot_id": "NONE",
            "best_composite_score": 0.0,
        }
    druggable = [e for e in evaluations if e.druggable_pocket_valid]
    best = max(evaluations, key=lambda e: e.composite_binding_score)
    return {
        "total_pockets": len(evaluations),
        "druggable_count": len(druggable),
        "druggable_fraction": round(len(druggable) / len(evaluations), 4),
        "mean_plddt": round(sum(e.pocket_mean_plddt for e in evaluations) / len(evaluations), 2),
        "best_uniprot_id": best.uniprot_id,
        "best_composite_score": round(best.composite_binding_score, 4),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="FoldEnergetics-xTB Structural & Ligand Triage")
    parser.add_argument("--uniprot", default="P00533")
    parser.add_argument("--min-plddt", type=float, default=75.0)
    args = parser.parse_args()

    triage = StructuralPocketTriage(min_plddt=args.min_plddt)
    plddts = [91.4, 88.7, 93.2, 86.5, 89.1, 90.8, 84.3]
    pae_vals = [2.1, 1.8, 3.2, 2.7, 4.1, 2.4, 3.8, 1.9]
    candidates = [
        CandidateLigand("CHEMBL939", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", 8.42, 1.15, 0.85, 3),
        CandidateLigand("CHEMBL553", "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", 7.95, 1.62, 1.10, 2),
        CandidateLigand("CHEMBL1946170", "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1O[C@H]1CCOC1", 8.85, 0.92, 0.72, 4),
    ]
    result = triage.score_pocket_and_ligands(args.uniprot, plddts, pae_vals, candidates)
    print(json.dumps(asdict(result), indent=2))


if __name__ == "__main__":
    main()
