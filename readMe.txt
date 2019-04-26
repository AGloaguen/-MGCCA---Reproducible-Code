Five folds in this directory and one file:

- matlab_toolboxes         : "N-way", "CMTF", "Poblano"and "Tensor" (please add them to your matlab path)
- pathdef.m                : please, once you have add 'matlab_toolboxes' to your path, save it to 'pathdef.m' and replace this file.
- RGCCA                    : the R package of RGCCA.
- Section_3.1.1_KKT        : code to generate Table 1 of Section 3.1.1
- Section_3.1.2_Simulation : code to generate Figure 1 and 2 of Section 3.1.2
- Section_3.2_EEG_data     : code to generate Figure 3 of Section 3.2

Three last folders also contain readMe.txt to explain how code works.

Required libraries for all R scripts are (please install them):

install.packages(c("devtools", "ROCR", "R.matlab", "grid", "gridExtra", "mvtnorm", "Matrix", "abind", "rjson", "tools", "ggplot2", "pracma", "pbmcapply", "multiway", "stringr", "CMA"))

WARNING : To install package CMA you have to execute the three following lines:

source("https://bioconductor.org/biocLite.R")
biocLite("CMA")

(taken from : https://bioconductor.riken.jp/packages/3.1/bioc/html/CMA.html)

Remark : EEG data are available in :

Section_3.2_EEG_data/tensor_EEG_data/tensor_EEG_4syllables_ISI_430_Subjdata_export_181121_2345_10_windows_normalisation_5_channels_rmved
