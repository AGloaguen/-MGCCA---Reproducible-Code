Step 1 : To run this code, open your terminal, go to this directory and execute "Rscript main.R JSON/data.json JSON/method.json"

Output : This will create a folder untitled "Results_EEG_data". Inside there is a folder untitle "Run_current_date_current_hour". Inside There is 2 folders : "figures" and "bootstrap_parameters_files" and 4 files : "cmd_Run_current_date_current_hour.txt", "reducers_parameter_file.Rdata", "cmtf_full_data_file.m" and "cmtf_full_data_parameter_file.mat".

Step 2 : Execute all lines of "cmd_Run_current_date_current_hour.txt" in your terminal in the same order.

Output : Generate a bunch of files. Results are present in "figures" where you can find a figure per method, per component and per mode.

Remark : Parameters are not set to generate the exact same figure as in the artcile in order for the simulation to be faster. To generate the exact same figure, set "bootstrap_samples" to "2000" (number of bootsrap samples) and "nb_of_bootstrap_group" to "50" (number of bootsrap sample treated per command line, not necessary).

Remark : On Windows, pbmclapply does not work and might cause some trouble. In that case, set "nb_cores" to "1" (number of cores to use for parallelization) in JSON/method.json.