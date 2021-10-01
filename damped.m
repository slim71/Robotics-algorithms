function damp_mat = damped(mat, lambda)
    n_row = size(mat, 1);
    mat_transpose = mat.';
    core = mat * mat_transpose + lambda^2 * eye(n_row);
    damp_mat = mat_transpose * core^(-1);
end
