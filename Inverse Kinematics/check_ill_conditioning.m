function [ill, num_cond] = check_ill_conditioning(mat)
% CHECK_ILL_CONDITIONING check if the supplied matrix is ill-conditioned.
%   The result is based on the conditioning number of the matrix, with a 
%   threshold chosen to be at 10^3.
%
% INPUT
% mat - matrix to check

    num_cond = cond(mat);
    ill = (num_cond >= 10^3);
end