from fold_energetics_xtb import CandidateLigand, StructuralPocketTriage


def test_pocket_and_ligand_triage() -> None:
    triage = StructuralPocketTriage(min_plddt=75.0)
    plddts = [91.4, 88.7, 93.2, 86.5, 89.1]
    pae_vals = [2.1, 1.8, 3.2, 2.7, 3.9]
    ligands = [
        CandidateLigand("CHEMBL_LOW", "CCO", 5.2, 4.8, 2.9, 1),
        CandidateLigand("CHEMBL_HIGH", "CC(=O)Nc1ccc(O)cc1", 8.7, 0.9, 0.75, 4),
    ]
    res = triage.score_pocket_and_ligands("P00533", plddts, pae_vals, ligands)
    assert res.druggable_pocket_valid is True
    assert res.top_ligand_id == "CHEMBL_HIGH"
    assert res.composite_binding_score > 10.0
