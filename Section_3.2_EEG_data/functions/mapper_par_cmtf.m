% parameters_file, output_file, file_number
run([pwd, '/../../pathdef.m'])
addpath(ans);
%% Parameters
load(parameters_file)
load(data_path)

modes{1}    = [1, 2, 3];
modes{2}    = [1, 4];
X1          = X1;
X2          = X2;
sz          = [size(X1) ,size(X2, 2)];
R           = nfac;

if isequal(init, 'svd')
    init = 'nvecs';
end

options              = ncg('defaults');
options.Display      ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = tol;
options.RelFuncTol   = tol;

if parfor_bool
    load(full_data_file)
    nb_boot = size(idx_boot, 1);
    Ystar  = zeros(sz(1), R, nb_boot);
    Astar  = zeros(sz(4), R, nb_boot);
    Bstar   = zeros(sz(2), R, nb_boot);
    Cstar   = zeros(sz(3), R, nb_boot);
    parfor (i = 1:nb_boot, nb_cores)
        idx             = idx_boot(i, :);
        Z               = par_load_tensors(X1, X2, modes, sz, idx);
        [Fac,G, out]    = cmtf_opt(Z, R, 'init', init , ...
            'alg_options',options); 
        Ystar(:, :, i) = Fac.U{1} * diag(sign(diag(corr(Fac.U{1}, Fac_full_data{1}))));
        Astar(:, :, i) = Fac.U{4} * diag(sign(diag(corr(Fac.U{4}, Fac_full_data{4}))));
        Bstar(:, :, i)  = Fac.U{2} * diag(sign(diag(corr(Fac.U{2}, Fac_full_data{2}))));
        Cstar(:, :, i)  = Fac.U{3} * diag(sign(diag(corr(Fac.U{3}, Fac_full_data{3}))));
    end
    save(output_file, 'Ystar', 'Astar', 'Bstar', 'Cstar', ...
        'file_number')
else
    idx           = 1:sz(1);
    Z             = par_load_tensors(X1, X2, modes, sz, idx);
    Fac_full_data = cmtf_opt(Z, R, 'init', init, 'alg_options', options); 
    save(output_file, 'Fac_full_data')
end
