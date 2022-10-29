function lbd = LieBracketsDistribution(dist1, dist2, vars)
% LIEBRACKETSDISTRIBUTION computes the Lie-Bracket between two
%   distributions.
%   This means applying the operator to each pair of vector fields in the
%   distributions.


    lbd = [];

    % Considering the number of columns for each distribution as the number
    % of vector fields in it
    for i = 1:size(dist1, 2)
        for j =1:size(dist2, 2)
            lb = LieBracket(dist1(:, i), dist2(:, j), vars);
            lbd = [lbd lb];
        end
    end
end
