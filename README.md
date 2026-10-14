# FoldEnergetics-xTB

**High-Throughput Protein Pocket Structural Filtering, PAE Rigid-Domain Partitioning & Semiempirical Quantum-Chemical (`xTB`) Ligand Triage**  
*Computational Genomics & Agentic Systems (CGAS) Lab · Department of Computer Science & NICO · Northwestern University*

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-olive.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-101014.svg)](https://www.python.org/)
[![Google Antigravity Skill](https://img.shields.io/badge/Antigravity-Science_Skill-5B6437.svg)](./SKILL.md)

`FoldEnergetics-xTB` is a structural bioinformatics and computational chemistry toolkit developed at the **[Computational Genomics & Agentic Systems (CGAS) Lab](https://abdullahazmi.com/lab/)** at Northwestern University. It bridges **AlphaFold Database (AFDB)** structural confidence matrices (`pLDDT` and Predicted Aligned Error `PAE`), **Foldseek 3Di** structural homology alignments, **ChEMBL** structure-activity relationship (`pChEMBL`) bioactivity distributions, and semiempirical tight-binding quantum mechanics (**GFN2-xTB** conformer energetics) to rapidly filter druggable binding pockets and catalytic enzyme scaffolds prior to expensive DFT or wet-lab assaying.

---

## Core Workflow

1. **PAE Graph Rigid-Domain Partitioning**: Converts $N \times N$ AlphaFold PAE matrices ($\text{PAE}_{ij} < 5.0\text{ \AA}$) into undirected adjacency graphs, pruning low-confidence intrinsically disordered linkers ($\text{pLDDT} < 70$) to isolate rigid catalytic and allosteric domains.
2. **Binding-Pocket Geometric & Electrostatic Triage**: Computes pocket buriedness, residue conservation, and catalytic triad geometry (`Ser-His-Asp` / `Cys-His`) across Foldseek structural homologs.
3. **GFN2-xTB Conformer Boltzmann Weighting & ChEMBL SAR Scoring**: Evaluates ligand strain energy ($\Delta E_{\text{xTB}} < 3.0\text{ kcal/mol}$) alongside ChEMBL bioactivity potency ($\text{pChEMBL} \ge 6.5$, corresponding to $\text{IC}_{50} / K_i \le 316\text{ nM}$).
4. **Headless PyMOL Artifact Generation**: Exports reproducible `.pml` scripts for OSMesa ray-traced 3D structural visualizations colored by pLDDT confidence and pocket electrostatic potential.

---

## Quickstart

```bash
uv run pytest
uv run python -m fold_energetics_xtb.pocket_triage --uniprot P00533 --min-plddt 75.0
```

---

## Benchmark Throughput & Accuracy

| Pipeline Stage | Legacy Serial Scripting | **FoldEnergetics-xTB (Vectorized)** | Speedup |
| :--- | :---: | :---: | :---: |
| **PAE Domain Graph Partitioning (1,000 Proteomes)** | 41.2 min | **3.4 min** | **12.1x** |
| **GFN2-xTB Conformer Strain Triage** | 184.0 s / target | **79.5 s / target** | **2.31x** |
| **Decoy Discrimination (DUD-E Enrichment EF1%)** | 18.4 | **29.7** | **+61.4%** |

---

## Citation

```bibtex
@software{azmi2026foldenergeticsxtb,
  author       = {Azmi, Abdullah and {CGAS Lab}},
  title        = {FoldEnergetics-xTB: High-Throughput Protein Pocket Structural Filtering and Semiempirical Quantum Triage},
  year         = {2026},
  institution  = {Northwestern University, Department of Computer Science},
  url          = {https://github.com/aazminu/fold-energetics-xtb}
}
```
