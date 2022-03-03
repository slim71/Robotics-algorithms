function dpinv = damped_pinv(mat, lambda)
    n = size(mat,1);
    core = mat * mat' + lambda^2 * eye(n);
    dpinv = mat' * inv(core);
end