function Cnew = RankOneUpdate(C, A, B)
% TODO

    % Supposed dimensions
    % matrix: m x n
    % u: m x k
    % v: k x n

    [m, n] = size(C);
    k = size(A, 2);

    Cupdated = C;

    for p = 1:k
        for j = 1:n
            for i = 1:m
                Cupdated(i, j) = Cupdated(i, j) + A(i, p) * B(p, j);
            end
        end
    end

    Cnew = Cupdated;

end