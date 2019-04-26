function [Z, INIT] = par_load_tensors(file_T1, file_T2, INIT_file, modes, sz)

    load(file_T1)
    load(file_T2)
    load(INIT_file)
    
    Z.modes        = modes;
    Z.size         = sz;
    Z.object{1}    = tensor(T1);
    Z.object{2}    = tensor(T2);
    INIT           = struct2cell(INIT);
    
end
