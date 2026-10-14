---
name: fold-energetics-xtb
description: Partition AlphaFold PAE matrices into rigid structural domains, score binding pocket confidence (pLDDT), and triage candidate ligands by GFN2-xTB conformer strain energy and ChEMBL pChEMBL potency.
---

# FoldEnergetics-xTB Agent Skill

Use this skill inside Google Antigravity alongside `google-deepmind/science-skills` (`alphafold_database_fetch_and_analyze`, `foldseek_structural_search`, `chembl_database`, `uniprot_database`, and `pymol`) for structural bioinformatics and small-molecule/enzyme design triage.

## Usage

```bash
uv run python -m fold_energetics_xtb.pocket_triage --uniprot P00533 --min-plddt 75.0
```
