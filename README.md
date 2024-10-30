# STged Tutorial and Example Guide

This tutorial provides a comprehensive guide on using STged with toy examples, covering both simulation and real data analyses.

## File Structure and Contents
The following folders contain all necessary code, data, and documentation to reproduce main and supplementary results in the STged study:

1. **Code Folder**:
- **Benchmark File**: Contains code for implementing benchmark methods used for comparison with STged.
- **Module File**: Contains code for STged-specific modules designed for analyzing real data.
- **Figures Code**: Includes scripts for generating figures from STged analyses.

2. **Simulation Folder**:
- Contains both **code** and **demo files** needed to reproduce results for main and supplementary STged simulation studies, enabling a full demonstration of STged's simulation-based performance.

3. **Demo Folder**:
- Contains **code** and **demo files** for reproducing STged results in real data analyses, specifically tailored to the examples used in the main and supplementary results sections.

## Tutorials
Detailed tutorials provide step-by-step guidance on running analyses on simulated and real data.

### Simulation Studies
These tutorials demonstrate STged’s functionality on simulated datasets, providing insight into its performance in controlled settings.
- **[seqFISH+ Data Analysis](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/Simulation/Demo-Simulation_result_FISH+.html)**: A demo showcasing STged applied to seqFISH+ data simulations.
- **[MERFISH Data Analysis](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/Simulation/Demo-Simulation_result_merfish.html)**: A tutorial on applying STged to MERFISH simulated data.

### Real Data Studies
Real data examples illustrate STged’s capabilities in practical applications, focusing on human PDAC tissue and mouse kidney tissue.
- **[mHVG/ctHVG Module Analysis](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/demo_files/demo_PDACA_STged_mHVG.html)**: Analyzing mHVG and ctHVG in human PDAC tissue, focusing on identifying highly variable genes linked to cell types.
- **[Cell Type Subpopulation Analysis](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/demo_files/demo_PDACA_STged_subpopulation.html)**: Identifying and analyzing cell subpopulations within human PDAC tissue.
- **[Cell-Cell Communication Analysis](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/demo_files/demo_PDAC_STged_CC.html)**: Analyzing cell-cell communication dynamics in human PDAC tissue.
- **[Gene Expression Program Module](https://htmlpreview.github.io/?https://github.com/TJJjiajuan/STged_example/blob/main/demo_files/demo_kidney_STged_downanalysis.html)**: Gene expression analysis of mouse kidney tissue, showcasing spatial gene expression programs.
