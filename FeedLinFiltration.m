function Gamma_fin = FeedLinFiltration(Delta, Gamma_0, vars)
% FEEDLINFILTRATION computes the Filtration of distributions Delta and Gamma_0.
%   This operation stops if, for some k, Gamma_k and Gamma_k+1 are
%   nonsingular at a point x and dim(Gamma_k(x)) = dim(Gamma_k+1(x)).
%   A distribution is said to be 'nonsingular' if, for all x in X, 
%   dim(Delta(x)) = n (system dimension).
%   Finally, since this computes a distribution, the function uses the
%   column range of a space to built the result.
    
    n = length(vars);

    Gamma_k = Gamma_0;
    Gamma_k1 = [];
    
    while (rank(Gamma_k) < n) && (rank(Gamma_k) ~= rank(Gamma_k1))
        Gamma_k1 = Gamma_k;

        lbd = LieBracketsDistribution(Delta, Gamma_k1, vars);

        % This batch of operations are equivalent to the sum of subspaces
        subspace = [Gamma_k1 lbd];
        Gamma_k = LinIndCols(subspace);
    end

    Gamma_fin = Gamma_k;
end
