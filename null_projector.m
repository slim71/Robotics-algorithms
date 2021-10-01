function n_proj = null_projector(mat)
    [~, columns] = size(mat);
    
    n_proj = eye(columns) - pinv(mat)*mat;
end