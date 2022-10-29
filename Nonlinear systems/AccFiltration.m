function minDelta = AccFiltration(delta, delta0, vars)
% ACCFILTRATION computes the Filtration of distributions to find the smallest
%   Delta-invariant distribution containing Delta0 (denoted <Delta, Delta0>).
%   This operation stops if, for some k, Delta_k and Delta_k+1 are
%   nonsingular at a point x and dim(Delta_k(x)) = dim(Delta_k+1(x)).
%   A distribution is said to be 'nonsingular' if, for all x in X, 
%   dim(Delta(x)) = n (system dimension).
%   Finally, since this computes a distribution, the function uses the
%   column range of a space to built the result.
    
    n = length(vars);

    delta_k = delta0;
    delta_k1 = [];
    
    while (rank(delta_k) < n) && (rank(delta_k) ~= rank(delta_k1))
        delta_k1 = delta_k;

        lbd = LieBracketsDistribution(delta_k1, delta, vars);

        % This batch of operations are equivalent to the sum of subspaces
        subspace = [delta_k1 lbd];
        delta_k = LinIndCols(subspace);
    end

    minDelta = delta_k;
end
