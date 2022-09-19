function obscod = ObservabilityCodistribution(codist, dist, vars)
% OBSERVABILITYSDISTRIBUTION computes the derivative of a covector field
%   along a vector field, between two distributions.
%   This means applying the operator to each (covector, vector) fields pair in 
%   the co/distributions.

    obscod = [];

    % Considering the number of columns/rows for the distribution/codistribution
    % as the number of vector/covector fields in it
    for i = 1:size(dist, 2)
        for j =1:size(codist, 1)
            covder = CovFieldDerivative(dist(:, i), codist(j, :), vars);
            obscod = [obscod; covder];
        end
    end
end
