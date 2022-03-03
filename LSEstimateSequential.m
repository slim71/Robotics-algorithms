function [Matpseudo, leftMat, rightMat] = LSEstimateSequential(matrix, row_num, lambda)

    if row_num >= size(matrix, 1) || row_num <= 0
        Matpseudo = Rank1UpdatePseudo(matrix);
        leftMat = Matpseudo;
        rightMat = [];
        return;
    end

    R = matrix(1:row_num, :);
    S = matrix(row_num+1:end, :);
    
    Rpseudo = pinv(R);%Rank1UpdatePseudo(R);%(R' * R)^(-1) * R';
    id_Rpseudo = eye(size(Rpseudo, 1));
    
    E = S * (id_Rpseudo - Rpseudo * R); 
%     E(E<1e-6)=0;
    
%     if min(size(E)) == 1 && svd
    Epseudo = damped_pinv(E, lambda);%Rank1UpdatePseudo(E);%(E' * E)^(-1) * E';
    id_E = eye(size(E, 1));
    id_Epseudo = eye(size(Epseudo, 1));
    Ecore = id_E - E * Epseudo;

    K = (id_E + Ecore * S * (Rpseudo * Rpseudo') * S' * Ecore)^(-1);

    T = Epseudo + (id_Epseudo - Epseudo * S) * (Rpseudo * Rpseudo') * ...
        S' * K * Ecore;
    leftMat = Rpseudo - T * S * Rpseudo;
    rightMat = T;
    Matpseudo = [leftMat, rightMat];
end