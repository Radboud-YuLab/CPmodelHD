# CP HD Model
---

High-dimensional CP (HD CP) is an extension of the original CP framework developed by the Yu group [^1]. It enables the analysis of bulk gene expression data across multiple genes simultaneously. This allows the model to capture the collective, stochastic interactions that drive cell-state transitions. 

Currently, HD CP models a single fate trajectory, making it suitable for linear or unidirectional processes. Future work will extend the framework to handle multifate or branching transitions in a heterogeneous cell population.

---

#### This framework includes:

Primary scripts: 
| Script | Description |
|--------|------------|
| CP_HD_main.m |Main script for running the HD CP model. |
| getTauHD.m | Calculate tau (parameter of a scaled progress ranging from 0 to 1) along the trajectory using distance in scaled time and expression space. |
| bootTimeSeriesHD.m  | Weighted bootstrapping to create evenly spaced timepoints and estimate gene expression using distance-weighted resampling. |
| getCpsHD.m | Fits the HD CP model to a given dataset and outputs CP coordinates. Time values are estimated using non-negative linear regression, and gene expressions are estimated using multiple linear regression. |
| getTopGenes.m | Extract top-expressed genes at each CP coordinate (excluding the first and last CP). |
| drawFigures3D.m | Generate 3D figures to visualize HD CP output. |








Reference:
[^1]: Chen CM, Yu R. 2025. A multi-step completion process model of cell plasticity. Briefings in Bioinformatics. 26(2). doi:https://doi.org/10.1093/bib/bbaf165. 
