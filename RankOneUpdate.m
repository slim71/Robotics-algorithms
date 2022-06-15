function Cnew = RankOneUpdate(C, A, B)
% RANKONEUPDATE computes the update of a matrix using the rank-one update
%   method.
%   The method is called "rank-one" because at each iteration only an
%   element of the matrix is updated (namely, a submatrix a | ranke(a)=1).
%   The rank-one update is used to efficiently compute the update of a
%   matrix, often in situations were this operation has to be done
%   iteratively.
%
%   Definition: Let M∈ℝ^(𝑚×𝑛) and column vectors 𝑢∈ℝ^𝑚 and 𝑣∈ℝ^𝑛. Then, 
%               the transformation M+𝑢𝑣^𝑇 is called a rank-one update to A.
%
%   In this instance, the method is not limited to only vector product: A
%   and B (u and v in the example) can indeed be entire matrices.
%   The method, not optimized here, computed the entire update and runs in
%   O(n^3).
%
% INPUT
%   C - original matrix (M in the example)
%   A - first product operand (u in the example)
%   B - second product operand (v in the example)


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