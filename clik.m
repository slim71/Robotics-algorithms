function qdot = clik(tt, q, J_sym, twist_sym)
    % initialization
    qdot = zeros(7,1);

    jac = J_sym(q);                 % 6x7
    p_jac = pinv(jac);              % 7x6
    twist_des = twist_sym(tt)';     % 6x1

    % actual function
    % qdot = J^(-1) * xdot_des
    for i=[1:size(qdot,1)]
        qdot(i,1) = p_jac(i,:) * twist_des;
    end
end