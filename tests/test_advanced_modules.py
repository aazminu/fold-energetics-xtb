import pytest
from fold_energetics_xtb import (
    ConformerState,
    Foldseek3DiHomologyScorer,
    GFN2xTBBoltzmannEvaluator,
    PAERigidDomainPartitioner,
    Residue3D,
)


def test_pae_rigid_domain_partitioning_two_domains_and_linker() -> None:
    # 8 residues: 0,1,2 in Domain 0; 3,4 disordered linker; 5,6,7 in Domain 1
    plddts = [92.0, 89.5, 91.0, 48.0, 52.0, 88.0, 90.5, 87.0]
    pae = [[12.0 for _ in range(8)] for _ in range(8)]
    for i in range(8):
        pae[i][i] = 0.5
    for i in (0, 1, 2):
        for j in (0, 1, 2):
            if i != j:
                pae[i][j] = 1.8
    for i in (5, 6, 7):
        for j in (5, 6, 7):
            if i != j:
                pae[i][j] = 2.1

    partitioner = PAERigidDomainPartitioner(min_plddt=70.0, max_pae_angstrom=4.5, min_domain_size=3)
    res = partitioner.partition(pae, plddts)

    assert len(res.domains) == 2
    assert res.domains[0].residue_indices == [0, 1, 2]
    assert res.domains[1].residue_indices == [5, 6, 7]
    assert res.disordered_linker_residues == [3, 4]
    assert res.interface_cross_pae[(0, 1)] == 12.0
    assert res.modularity_score > 0.35


def test_pae_partition_edge_cases() -> None:
    partitioner = PAERigidDomainPartitioner()
    with pytest.raises(ValueError):
        partitioner.partition([], [])
    with pytest.raises(ValueError):
        partitioner.partition([[1.0, 2.0]], [90.0, 90.0])
    with pytest.raises(ValueError):
        partitioner.partition([[-1.0]], [90.0])

    # All disordered
    res = partitioner.partition([[1.0, 1.0], [1.0, 1.0]], [40.0, 50.0])
    assert res.domains == []
    assert res.disordered_linker_residues == [0, 1]


def test_foldseek_3di_and_catalytic_triad_detection() -> None:
    scorer = Foldseek3DiHomologyScorer()
    pocket = [
        Residue3D(195, "SER", (0.0, 0.0, 0.0), (0.5, 0.5, 0.5)),
        Residue3D(57, "HIS", (3.2, 0.0, 0.0), (3.4, 0.4, 0.2)),
        Residue3D(102, "ASP", (6.1, 0.5, 0.0), (6.3, 0.6, 0.3)),
        Residue3D(214, "GLY", (9.5, 1.2, 0.5), (9.5, 1.2, 0.5)),
    ]
    seq_3di = scorer.encode_3di_sequence(pocket)
    assert len(seq_3di) == 4

    raw_score, identity = scorer.align_3di(seq_3di, seq_3di)
    assert identity == 1.0
    assert raw_score > 10.0

    # Shifted sequence with a single N-terminal insertion must recover high aligned identity via NW traceback
    gapped_score, gapped_identity = scorer.align_3di("ACDEFGH", "CACDEFGH")
    assert gapped_score == 22.5
    assert gapped_identity == 0.875

    with pytest.raises(ValueError, match="valid 3Di tokens"):
        scorer.align_3di("ACDEF", "ACXEF")

    # Translation invariance of dRMSD
    translated = [
        Residue3D(
            r.res_idx,
            r.res_name,
            (r.ca_xyz[0] + 15.0, r.ca_xyz[1] - 8.0, r.ca_xyz[2] + 3.0),
            (r.cb_xyz[0] + 15.0, r.cb_xyz[1] - 8.0, r.cb_xyz[2] + 3.0),
        )
        for r in pocket
    ]
    assert scorer.compute_distance_rmsd(pocket, translated) == 0.0

    triads = scorer.detect_catalytic_triads(pocket)
    assert len(triads) == 1
    assert triads[0].residues == (195, 57, 102)
    assert "Ser-His-Asp" in triads[0].triad_type


def test_gfn2_xtb_boltzmann_ensemble_and_pymol_script() -> None:
    evaluator = GFN2xTBBoltzmannEvaluator(temperature_kelvin=298.15)
    conformers = [
        ConformerState("CONF_MIN", -42.500000, -6.20, 0.65, 4),
        ConformerState("CONF_STRAINED", -42.497500, -5.90, 1.25, 2),
        ConformerState("CONF_HIGH", -42.492000, -5.10, 2.40, 1),
    ]
    res = evaluator.evaluate_ensemble(conformers, bioactive_conformer_id="CONF_STRAINED")
    assert res.lowest_energy_conformer == "CONF_MIN"
    assert abs(sum(res.boltzmann_weights.values()) - 1.0) < 1e-4
    assert res.boltzmann_weights["CONF_MIN"] > res.boltzmann_weights["CONF_STRAINED"]
    assert res.free_energy_penalty_kcal_mol > 0.0

    pml = evaluator.generate_pymol_pocket_script("P00533", [57, 102, 195], "CHEMBL939")
    assert "select catalytic_pocket, resi 57+102+195" in pml
    assert "AF-P00533-F1-model_v4" in pml


def test_cli_benchmark_markdown_json_and_reproducibility(
    capsys: pytest.CaptureFixture[str], tmp_path: pytest.TempPathFactory
) -> None:
    import json
    import subprocess
    import sys
    from pathlib import Path
    from fold_energetics_xtb.__main__ import format_markdown_table, main, run_benchmark

    rows_a = run_benchmark(seed=42)
    rows_b = run_benchmark(seed=42)
    assert rows_a == rows_b
    assert len(rows_a) == 4

    md = format_markdown_table(rows_a)
    assert "| benchmark | target | primary_metric | score |" in md
    assert "pocket_ligand_triage" in md
    assert "gfn2_xtb_boltzmann_ensemble" in md

    out_file = Path(str(tmp_path)) / "bench.json"
    rc = main(["--format", "json", "--seed", "42", "--output", str(out_file)])
    assert rc == 0
    captured = capsys.readouterr().out
    payload = json.loads(captured)
    assert payload["seed"] == 42
    assert payload["benchmarks"] == rows_a
    assert json.loads(out_file.read_text(encoding="utf-8")) == payload

    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, str(repo_root / "scripts" / "benchmark.py"), "--format", "markdown"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "pae_rigid_domain_partition" in proc.stdout

    with pytest.raises(SystemExit):
        main(["--format", "invalid_fmt"])

