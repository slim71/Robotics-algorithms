function lin_ind_cols = LinIndCols(matrix)
% LININDCOLS extract the linerarly independent columns of the supplied
%   matrix.

    retmat = matrix(:,1);
    
    % Add one column at a time and see if rank increases
    for i = 2:size(matrix, 2)
       if rank([retmat matrix(:, i)]) > rank(retmat)
           retmat = [retmat matrix(:, i)];
       end
    end
    
    lin_ind_cols = retmat;
    
end
