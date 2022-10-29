function minGamma = ObsFiltration(delta, gamma0, vars)
% OBSFILTRATION computes the Filtration of a codistribution Gamma0 and a
%   distribution delta to find the smallest Delta-invariant codistribution
%   which contains Gamma0 (denoted <Delta, Gamma0>).
%   This operation stops if, for some k, Gamma_k and Gamma_k+1 are
%   nonsingular at a point x and dim(Gamma_k(x)) = dim(Gamma_k+1(x)).
%   A distribution is said to be 'nonsingular' if, for all x in X, 
%   dim(Delta(x)) = n (system dimension).
%   Finally, since this computes a codistribution, the function uses the
%   row range of a space to built the result.

    n = length(vars);

    gamma_k = gamma0;
    gamma_k1 = [];

    % A distribution is said involutive if the Lie-bracket between any of its 
    % vector fields remains in the distribution
    while (rank(gamma_k) < n) && (rank(gamma_k) ~= rank(gamma_k1))
        gamma_k1 = gamma_k;

        lbd = ObservabilityCodistribution(gamma_k1, delta, vars);

        % This batch of operations are equivalent to the sum of subspaces
        subspace = [gamma_k1; lbd];
        gamma_k = LinIndRows(subspace);
    end

    minGamma = gamma_k;
end
