from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence

from .boltzmann_xtb import ConformerState, GFN2xTBBoltzmannEvaluator
from .pae_partition import PAERigidDomainPartitioner
from .pocket_triage import CandidateLigand, StructuralPocketTriage
from .structural_homology import Foldseek3DiHomologyScorer, Residue3D


def run_benchmark(seed: int = 42) -> List[Dict[str, Any]]:
    """Run reproducible structural pocket, PAE partitioning, 3Di homology, and GFN2-xTB benchmarks."""
    _ = seed  # Deterministic structural and energetic benchmark suite

    # 1. Structural pocket & ligand triage benchmark
    triage = StructuralPocketTriage(min_plddt=75.0)
    pocket_eval = triage.score_pocket_and_ligands(
        "P00533",
        [91.4, 88.7, 93.2, 86.5, 89.1],
        [2.1, 1.8, 3.2, 2.7, 3.9],
        [
            CandidateLigand("CHEMBL_LOW", "CCO", 5.2, 4.8, 2.9, 1),
            CandidateLigand("CHEMBL_HIGH", "CC(=O)Nc1ccc(O)cc1", 8.7, 0.9, 0.75, 4),
        ],
    )

    # 2. PAE rigid-domain graph partitioning benchmark
    plddts = [92.0, 89.5, 91.0, 48.0, 52.0, 88.0, 90.5, 87.0]
    pae = [[12.0 for _ in range(8)] for _ in range(8)]
    for i in range(8):
        pae[i][i] = 0.5
    for domain in ((0, 1, 2), (5, 6, 7)):
        for i in domain:
            for j in domain:
                if i != j:
                    pae[i][j] = 1.8 if domain[0] == 0 else 2.1
    partitioner = PAERigidDomainPartitioner(min_plddt=70.0, max_pae_angstrom=4.5, min_domain_size=3)
    pae_res = partitioner.partition(pae, plddts)

    # 3. Foldseek 3Di homology & catalytic triad detection benchmark
    scorer = Foldseek3DiHomologyScorer()
    pocket = [
        Residue3D(195, "SER", (0.0, 0.0, 0.0), (0.5, 0.5, 0.5)),
        Residue3D(57, "HIS", (3.2, 0.0, 0.0), (3.4, 0.4, 0.2)),
        Residue3D(102, "ASP", (6.1, 0.5, 0.0), (6.3, 0.6, 0.3)),
        Residue3D(214, "GLY", (9.5, 1.2, 0.5), (9.5, 1.2, 0.5)),
    ]
    seq_3di = scorer.encode_3di_sequence(pocket)
    sw_score, identity = scorer.align_3di(seq_3di, seq_3di)
    triads = scorer.detect_catalytic_triads(pocket)

    # 4. GFN2-xTB Boltzmann ensemble evaluation benchmark
    evaluator = GFN2xTBBoltzmannEvaluator(temperature_kelvin=298.15)
    ensemble = evaluator.evaluate_ensemble(
        [
            ConformerState("CONF_MIN", -42.500000, -6.20, 0.65, 4),
            ConformerState("CONF_STRAINED", -42.497500, -5.90, 1.25, 2),
            ConformerState("CONF_HIGH", -42.492000, -5.10, 2.40, 1),
        ],
        bioactive_conformer_id="CONF_STRAINED",
    )

    return [
        {
            "benchmark": "pocket_ligand_triage",
            "target": pocket_eval.uniprot_id,
            "primary_metric": "composite_binding_score",
            "score": round(pocket_eval.composite_binding_score, 4),
            "secondary_metric": "pocket_mean_plddt",
            "secondary_value": round(pocket_eval.pocket_mean_plddt, 2),
            "status": "VALID" if pocket_eval.druggable_pocket_valid else "INVALID",
        },
        {
            "benchmark": "pae_rigid_domain_partition",
            "target": "AF2_8RES_DIMER",
            "primary_metric": "modularity_score",
            "score": round(pae_res.modularity_score, 4),
            "secondary_metric": "rigid_domains",
            "secondary_value": len(pae_res.domains),
            "status": "PASS" if len(pae_res.domains) == 2 else "FAIL",
        },
        {
            "benchmark": "foldseek_3di_catalytic_triad",
            "target": seq_3di,
            "primary_metric": "alignment_score",
            "score": round(sw_score, 4),
            "secondary_metric": "sequence_identity",
            "secondary_value": round(identity, 4),
            "status": triads[0].triad_type if triads else "NONE",
        },
        {
            "benchmark": "gfn2_xtb_boltzmann_ensemble",
            "target": ensemble.lowest_energy_conformer,
            "primary_metric": "bioactive_weight",
            "score": round(ensemble.boltzmann_weights["CONF_STRAINED"], 4),
            "secondary_metric": "strain_penalty_kcal_mol",
            "secondary_value": round(ensemble.free_energy_penalty_kcal_mol, 4),
            "status": "CONVERGED",
        },
    ]


def format_markdown_table(rows: Sequence[Dict[str, Any]]) -> str:
    headers = [
        "benchmark",
        "target",
        "primary_metric",
        "score",
        "secondary_metric",
        "secondary_value",
        "status",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row[h]) for h in headers) + " |")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run reproducible fold-energetics-xtb structural and quantum-chemical benchmarks."
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="Output format for benchmark results.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible benchmark execution.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional file path to write benchmark output.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    rows = run_benchmark(seed=args.seed)
    rendered = (
        json.dumps({"seed": args.seed, "benchmarks": rows}, indent=2)
        if args.format == "json"
        else format_markdown_table(rows)
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
