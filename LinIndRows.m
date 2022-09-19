function lin_ind_rows = LinIndRows(matrix)
% LININDCOLS extract the linerarly independent row of the supplied
%   matrix.

    retmat = matrix(1, :);
    
    % Add one row at a time and see if rank increases
    for i = 2:size(matrix, 1)
       if rank([retmat; matrix(i, :)]) > rank(retmat)
           retmat = [retmat; matrix(i, :)];
       end
    end
    
    lin_ind_rows = retmat;
    
end
