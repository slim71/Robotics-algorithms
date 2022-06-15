function ind = check_independence(a, b)
% CHECK_INDEPENDENCE combines the provided arrays to check if they are
%   independent one another.
%   The check system is based on the rank of a virtual matrix built
%   concatenating the provided arrays.
%
% INPUT
%   a - first array
%   b - second array

    together = [a, b];
    ind = ~(rank(together) < min(size(together)));
end