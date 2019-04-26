#################################
###    Folders description    ###
#################################

There are 4 folds in this directory:

- RGCCA                    : the R package of RGCCA.
- Section_3.1.1_KKT        : code to generate Table 1 of Section 3.1.1
- Section_3.1.2_Simulation : code to generate Figure 1 and 2 of Section 3.1.2
- Section_3.2_EEG_data     : code to generate Figure 3 of Section 3.2

Three last folders also contain readMe.txt to explain how code works.

#################################
###     Matlab libraries      ###
#################################

Please download the following 4 matlab toolboxes (with the right version):

- v3.2 of the "N-way" toolbox   : https://www.mathworks.com/matlabcentral/fileexchange/1088-the-n-way-toolbox?s_tid=mwa_osa_a
- v1.1 of the "CMTF" toolbox    : http://www.models.life.ku.dk/joda/CMTF_Toolbox
- v1.1 of the "Poblano" toolbox : https://software.sandia.gov/trac/poblano/downloader/download/file/15/poblano_toolbox_1.1.zip
- v2.6 of the "Tensor" toolbox  : https://www.sandia.gov/~tgkolda/TensorToolbox/reg-2.6.html

Once you have done so, please add these 4 toolboxes to your path and save it to 'pathdef.m' (replace the current 'pathdef.m' file).

#################################
###       R libraries         ###
#################################

Required libraries for all R scripts are (please install them):

install.packages(c("devtools", "ROCR", "R.matlab", "grid", "gridExtra", "mvtnorm", "Matrix", "abind", "rjson", "tools", "ggplot2", "pracma", "pbmcapply", "multiway", "stringr", "CMA"))

WARNING : To install package CMA you have to execute the three following lines:

source("https://bioconductor.org/biocLite.R")
biocLite("CMA")

(taken from : https://bioconductor.riken.jp/packages/3.1/bioc/html/CMA.html)

#################################
###         Various           ###
#################################

Remark : EEG data are available in :

Section_3.2_EEG_data/tensor_EEG_data/tensor_EEG_4syllables_ISI_430_Subjdata_export_181121_2345_10_windows_normalisation_5_channels_rmved
