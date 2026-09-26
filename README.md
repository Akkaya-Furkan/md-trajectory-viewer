# MD trajectory viewer

**Question.** A browser tool to load, play and analyse molecular dynamics trajectories.

## Data
Topology + trajectory files (PDB/GRO + XTC/DCD); demo data from the egfr-mmgbsa runs.

## Approach
1. Streamlit + MDAnalysis + py3Dmol
2. Frame slider with 3D view; live RMSD and radius of gyration (Plotly)
3. Per-residue RMSF and MDAnalysis selection syntax
4. Striding/subsampling and upload limits for large files; pytest for analysis functions

## Deliverables
Deployed app with a demo trajectory and a GIF in this README.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
