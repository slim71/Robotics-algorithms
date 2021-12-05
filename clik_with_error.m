function qdot = clik_with_error(tt, q, J_sym, twist_sym, err_sym, gains)
    % initialization
    qdot = zeros(7,1);

    jac = J_sym(q);                 % 6x7
    p_jac = pinv(jac);              % 7x6
    twist_des = twist_sym(tt)';     % 6x1
    error = err_sym(tt, q)';        % 6x1
    K_matrix = gains*eye(6);        % 6x6

    % actual function
    % qdot = J^(-1) * (xdot_des + K * e)
    for i=[1:size(qdot,1)]
        qdot(i,1) = p_jac(i,:) * (twist_des + K_matrix * error);
    end
end