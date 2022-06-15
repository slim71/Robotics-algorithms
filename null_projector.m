function n_proj = null_projector(mat)
% NULL_PROJ simply computes the Null Projector of a matrix using a
%   standard formula.
%
% INPUT
%   mat - which matrix to compute the projector of

    [~, columns] = size(mat);
    
    n_proj = eye(columns) - pinv(mat)*mat;
end