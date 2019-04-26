Step 1 : To run this simulation two possibilities :

	1) Using command lines
		i) Go to this directory
		ii) Execute in your Terminal "Rscript main.R"
	2) Open R or Rstudio (define working directory at location of this file) and run main.R

Output : This will create a folder untitled "simu_current_date_and_hour" in which 2 files are present : "commands.txt" and "mapper_cmtf_parameter_file.mat".

Step 2 : Execute all lines of "commands.txt" in your terminal in the same order.

Output : Generate a bunch of folders. Results are present in "figures" where "All_Methods_comp_separated_NoAxisTitle_DB.pdf" is Figure 1 of Section 3.1.2 and "fig_Weights_OneTensor_OneMode_SNR_0.5.pdf" is Figure 2.

Remark : Parameters are not set to generate the exact same figure as in the article in order for the simulation to be faster. To generate the exact same figure, set "n_simu" to "100" (number of data set generated), "n_random_start" to "10" (number of random start for the first component, used for MGCCA and RGCCA), "nstart_next_comps" to "10" (number of random start for the second component, used for MGCCA and RGCCA) and "SNR" to "c(0.1, 0.2, 0.3, 0.5, 0.7, 1)" (number of SNR evaluated).
