
run([pwd, '/../pathdef.m'])
addpath(ans);

%% Parameters
modes{1}    = [1, 2, 3];
modes{2}    = [1, 4, 5];
sz          = [50, 30, 30, 30, 30];

pathtoReadfiles      = strcat(path, '/SNR_', num2str(SNR), '/matfiles/');
fileList             = dir(strcat(pathtoReadfiles, 'Simu_Fold_*_T1.mat'));
nb_files             = size(fileList, 1);

files_to_rm     = dir(strcat(pathtoReadfiles, '../CMTF/cmtf.best_start_Fold_*.mat'));
idx_files_to_rm = [];

for i = 1:length(files_to_rm)
    idx_files_to_rm = [idx_files_to_rm, sscanf(files_to_rm(i, 1).name,'cmtf.best_start_Fold_%d.mat')];
end

files_to_treat = setdiff(1:nb_files, idx_files_to_rm);

parfor (i = files_to_treat, nb_cores)
    options              = ncg('defaults');
    options.Display      ='final';
    options.MaxFuncEvals = 100000;
    options.MaxIters     = 10000;
    options.StopTol      = tol;
    options.RelFuncTol   = tol;
    
    files     = dir(strcat(pathtoReadfiles, 'Simu_Fold_', num2str(i), '_T*.mat'));
    INIT_file = dir(strcat(pathtoReadfiles, 'INIT_simu_Fold_', num2str(i), '.mat'));
    
    file_T1   = strcat(pathtoReadfiles, files(1, 1).name);
    file_T2   = strcat(pathtoReadfiles, files(2, 1).name);
    INIT_file = strcat(pathtoReadfiles, INIT_file.name);
    [Z, INIT] = par_load_tensors(file_T1, file_T2, INIT_file, modes, sz);
    INIT      = cellfun(@(x) x(:, 1:R),INIT, 'UniformOutput',false);
    
    [best_Fac,best_G, best_out]   = cmtf_opt(Z, R, 'init', INIT ,'alg_options',options); 
    best_start                    = 'svd';
        
    % fit CMTF-OPT
    for j = 2:nstart
        [Fac,G,out]   = cmtf_opt(Z, R, 'init', 'random','alg_options',options); 
        if out.Fit > best_out.Fit
            best_Fac   = Fac;
            best_G     = G; 
            best_out   = out;
            best_start = 'random';
        end
    end
    
    Fac = best_Fac;
    
    parsave(strcat(path, '/SNR_', num2str(SNR), '/CMTF/cmtf.results_nbcomp_', num2str(R), '_Fold_', num2str(i), '.mat'), Fac)
    parsave(strcat(path, '/SNR_', num2str(SNR), '/CMTF/cmtf.best_start_nbcomp_', num2str(R), '_Fold_', num2str(i), '.mat'), best_start)

end
