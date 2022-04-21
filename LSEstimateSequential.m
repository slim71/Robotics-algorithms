function [Matpseudo, leftMat, rightMat] = LSEstimateSequential(matrix, row_num, lambda)
% TODO help

    % In case the subdivision is not valid
    if row_num >= size(matrix, 1) || row_num <= 0
        Matpseudo = blockMatPseudo(matrix);
        leftMat = Matpseudo;
        rightMat = [];
        return;
    end

    % Initialization
    R = matrix(1:row_num, :);
    S = matrix(row_num+1:end, :);

    % Round R for accuracy and to improve performance
    R = round(R, 5);

    % If R is row/column matrix, cannot check conditioning number properly
    if min(size(R)) == 1 || cond(R) > 1e3
        Rpseudo = damped_pinv(R, lambda);
    else
        % Tried but always too slow for practical use
        % Rpseudo = blockMatPseudo(R);
        Rpseudo = pinv(R);
    end

    id_Rpseudo = eye(size(Rpseudo, 1));
    E = S * (id_Rpseudo - Rpseudo * R); 
    
    % row/column matrix, cannot check conditioning number
    if min(size(E)) == 1 || cond(E) > 1e3
        Epseudo = damped_pinv(E, lambda);
    else
        Epseudo = blockMatPseudo(E, lambda);
        % As comparison, you can check numerical solution differences using
        % Epseudo = pinv(E);
    end

    id_E = eye(size(E, 1));
    id_Epseudo = eye(size(Epseudo, 1));
    Ecore = id_E - E * Epseudo;

    K = (id_E + Ecore * S * (Rpseudo * Rpseudo') * S' * Ecore)^(-1);

    T = Epseudo + (id_Epseudo - Epseudo * S) * (Rpseudo * Rpseudo') * ...
        S' * K * Ecore;

    % Build output matrices
    leftMat = Rpseudo - T * S * Rpseudo;
    rightMat = T;
    Matpseudo = [leftMat, rightMat];
end