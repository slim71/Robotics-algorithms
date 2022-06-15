function qdot = rp(tt, qq, sym_Js, sym_twistes, lambda, with_rou)
% RP solves inverse kinematics with the Reverse Priority method.
%   The Reverse Priority method was depicted by F. Flacco and A. De Luca in
%   their "A reverse priority approach to multi-task control of redundant
%   robots".
%
% INPUT
%   tt          - (int) time istant
%   qq          - (int) joint variables array
%   sym_Js      - (handle) Jacobian function handle
%   sym_twistes - (handle) array of twist function handles
%   lambda      - (int) dampening factor
%   with_rou    - (boolean) flag to use the Rank-one update

    % Initialization
    qdot_lp1 = zeros(max(size(qq)),1);

    qdot_k = qdot_lp1;
    JRA_k = [];

    for k = max(size(sym_twistes)):-1:1
        Jk = sym_Js{k}; % Jacobian function handle
        twistk = sym_twistes{k}; % Task function handle

        % Computed Jacobian
        Jk_qq = Jk(qq);
        % Computing Augmented Jacobian
        JRA_k = [Jk_qq; JRA_k];

        [~, Tk, ~] = LSEstimateSequential(JRA_k, size(Jk_qq, 1), lambda, with_rou);


        % Solution for k=l...1:
        % q'_k = q'_(k+1) + T_k * inv(J_k*T_k) * (task_k - J_k*q'_(k+1))
        qdot_k = qdot_k + Tk * pinv(Jk_qq*Tk) * (twistk(tt,qq)' - Jk_qq * qdot_k);
    end

    % Arrange solution vector for ODE solvers
    qdot = zeros(max(size(qq)),1);
    for k = 1:max(size(qq))
        qdot(k, 1) = qdot_k(k);
    end
end