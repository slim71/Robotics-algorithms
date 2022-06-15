function qdot = clik(tt, qq, J_sym, twist_sym)
% CLIK solves inverse kinematics with the classic Closed Loop method.
%   A good reference to study this method can be "Robotics: Modelling,
%   Planning and Control" by L. Sciavicco, L. Villani, G. Oriolo, B.
%   Siciliano.
%
% INPUT
%   tt          - (int) time istant
%   qq          - (int) joint variables array
%   J_sym       - (handle) Jacobian function handle
%   twist_sym   - (handle) twist function handle

    % initialization
    qdot = zeros(7,1);

    jac = J_sym(qq);                % 6x7
    p_jac = pinv(jac);              % 7x6
    twist_des = twist_sym(tt)';     % 6x1

    % actual function
    % qdot = J^(-1) * xdot_des
    for i= 1:size(qdot,1)
        qdot(i,1) = p_jac(i,:) * twist_des;
    end
end