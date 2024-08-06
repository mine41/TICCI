# Trajectory Inference with Cell–Cell Interactions (TICCI): intercellular communication improves the accuracy of trajectory inference methods

#### Author: Yifeng Fu

### This R and Python script implements the TICCI algorithm, which leveraging intercellular communication information examining cell trajectories in scRNA-seq data.

### Environment Configuration:

#### Our script is mainly based on R (version 4.32) and Python (version 3.7.16) with several packages as follows:

- Python: numpy, pandas, matplotlab, scanpy, anndata, sklearn, etc.
- R: cluster, Seurat, dplyr, ggplot2, CellChat v1, ggalluvial, svglite, scran, etc.
- For packages required but not list above, please refer to corresponding website for installation.

### Directory Overview:

- Data: Dataset
- Figures: Result images
- Output: Each output result
- Python tools: Encapsulated Python functions
- R_tools: Encapsulated R functions
- Rdata: Save records of variables in each stage of running R

### Scripts in the SLICE package

1. AT2_1.R: Preprocessing data and generating CCI information on AT2 data
2. AT2_2.ipynb: Reconstructing trajectories and Analyzing pseudotime on AT2 data
3. HSMM_1.R: Preprocessing data and generating CCI information on HSMM data
4. HSMM_2.ipynb: Reconstructing trajectories and Analyzing pseudotime on HSMM data

### TICCI is based on following datasets which are downloaded from Guo M, Bao EL, Wagner M, Whitsett JA, Xu Y. 2016:

Dataset can be downloaded from: http://research.cchmc.org/pbge/slice.html

1. GSE52583.Rda: The R data object file for the SLICE demonstration contains scRNA-seq data from Treutlein et al., 2014. It includes alveolar type 2 (AT2) cells from E14.5, E16.5, E18.5, and adult mouse lungs, as well as three populations of epithelial cells at E18.5 (EPI). Expression profiles with identical gene symbol annotations were averaged.
2. hs_km.Rda: This R data object file contains a pre-compiled pairwise Kappa similarity matrix for human genes. The similarity is based on GO_BP_FAT annotations, which were downloaded from DAVID on December 13, 2015.
3. HSMM.Rda: This R data object file contains scRNA-seq data from differentiating human skeletal muscle myoblasts (HSMM), as presented in Trapnell et al., 2014. Expression profiles with identical gene symbol annotations were averaged.
4. mm_km.Rda: This R data object file contains a pre-compiled pairwise Kappa similarity matrix for mouse genes. The similarity is based on GO_BP_FAT annotations, which were downloaded from DAVID on December 13, 2015.

### Implementation:

1. Download the datasets from the website: http://research.cchmc.org/pbge/slice.html.
2. Put the four datasets into the folder data: GSE52583.Rda, hs_km.Rda, HSMM.Rda, mm_km.Rda.
3. Open R GUI and install the necessary dependencies.
4. Run the R scripts.
   - Input: RDS raw data.
   - Output: Preprocessed data_csv, state label, single cell entropy, kmeans label, cellchat_net result
5. Run the root directory Python scrpits.
   - Input: preprocessed data_csv, state label, single cell entropy, kmeans label, cellchat_net result
   - Output: Final result h5ad, result images
