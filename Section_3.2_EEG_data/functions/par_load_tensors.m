function Z = par_load_tensors(X1, X2, modes, sz, idx)

    Z.modes     = modes;
    Z.size      = sz;
    Z.object{1} = tensor(centering(X1(idx, :, :), 'c'));
    Z.object{2} = tensor(double(X2(idx, :) - mean(X2(idx, :))));
        
end
