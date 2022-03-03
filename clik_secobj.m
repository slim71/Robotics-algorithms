function qdot = clik_secobj(tt, q, J_sym, twist_sym, err_sym, gains, sec_obj)
    % initialization
    qdot = zeros(7,1);

    jac = J_sym(q);                 % 6x7
    p_jac = pinv(jac);              % 7x6
    twist_des = twist_sym(tt)';     % 6x1
    error = err_sym(tt, q)';        % 6x1
    K_matrix = gains*eye(6);        % 6x6

    proj = eye(7) - p_jac*jac;      % 7x7
%     k0 = eye(7);%K_matrix; % for now
%     sec_obj = sqrt(det(jac*jac'));
%     
%     q0dot = k0 * gradient(sec_obj,q);
    q0dot = sec_obj(q);             % 7x1

    % actual function
    % qdot = J^(-1) * (xdot_des + K * e)
    for i=[1:size(qdot,1)]
        qdot(i,1) = p_jac(i,:) * (twist_des + K_matrix * error) + proj(i,:) * q0dot;
    end
end