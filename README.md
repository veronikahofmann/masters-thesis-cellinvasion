# masters-thesis-cellinvasion

This project was developed in line with my Master's Thesis about cell invasion models. It is an extension of the model from René F. M. van Oers et al., "Mechanical Cell-Matrix Feedback Explains Pairwise and Collective Endothelial Cell Behavior In Vitro," PLOS Computational Biology vol. 10, no. 8, 2014. https://doi.org/10.1371/journal.pcbi.1003774.

To use it, you need Python 3.9 or newer. It possibly also works with older versions, but I did not test it with them.

The simulation can be started by executing the file main.py. It will use the parameters specified in parameters.py and produce some animated mp4-plots as well as a line plot of the cell- and ECM-concentration. The plots will be saved as pngs in a folder named "PLOTS", the results will be saved as .out-files in a folder called "RESULTS".

Additionally to the CPM-FEM model, the modules colsons_eq.py and crossleys_eq.py provide implementations and plotting functions of Chloé Colson's parameterized PDE model (from C. Colson et al., "Travelling-wave analysis of a model of tumour invasion with degenerate, cross-dependent diffusion," Proceedings. Mathematical, physical, and engineering sciences, vol. 477, no. 2256, p. 20210593, 2021. https://doi.org/10.1098/rspa.2021.0593) and Rebecca Crossley's dimensionless PDE model (from R. M. Crossley et al., "Travelling waves in a coarse-grained model of volume-filling cell invasion: Simulations and comparisons," Studies in Applied Mathematics, vol. 151, pp. 1471-1497, 2023. https://doi.org/10.1111/sapm.12635), respectively.
