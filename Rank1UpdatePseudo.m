function pseudo = Rank1UpdatePseudo(matrix)

    % A contains an incrementing number of columns of the supplied matrix,
    % becoming a submatrix of it
    A = matrix(:, 1); % A_1, also a_1 here

    n_col = size(matrix, 2);

    % Pseudoinverse of the submatrix: (A_1)^†
    if all(A == 0) % useless?
        Apseudo = 0; % zeros(1, size(matrix,1));
    else
        Apseudo = (A' * A)^(-1) * A';
    end

    for k = [2:n_col]
        % k-th column of the matrix to invert
        ak = matrix(:, k); % a_k

        % d_k = A_(k−1)^† * a_k
        dk = Apseudo * ak;

        % c_k = a_k − A_(k−1) * d_k
        ck = ak - A * dk;

        % c_k=0,  b_k = (1 + d_k^T * d_k)^(−1) * d_k^T * A_(k−1)^†
        if all(ck == 0)
            bk = (1 + dk' * dk)^(-1) * dk' * Apseudo;
        else % c_k≠0,  b_k = c_k^†
            bk = (ck' * ck)^(-1) * ck';
        end

	    % A_k^†= ( (A_(k−1)^† − d_k * b_k) ¦  b_k )
        Apseudo = [Apseudo - dk * bk; bk];

        % Increase submatrix A
        A = matrix(:, 1:k);

    end

    pseudo = Apseudo;

end