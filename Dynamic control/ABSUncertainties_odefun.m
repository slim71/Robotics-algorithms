function unc_pars_dot = ABSUncertainties_odefun(t, y, q_array, des_array, lambda, R)
    
    n = size(q_array, 2)/2; % or size(des_array, 2)

    qs = q_array(1:n);
    dqs = q_array(n+1:2*n);

    % Simply any desired configuration derived from Cinematic Inversion
    q_des = des_array(1:n)';
    dq_des = des_array(n+1:2*n)';
    ddq_des = des_array(2*n+1:end)';

    % Errors
    e = q_des - qs;
    de = dq_des - dqs;
    s = de' + lambda*e';

    % Regressor
    Y = DynamicRegressor(qs, dqs, dq_des, ddq_des);
    
    % Update law (u_pi) of estimated parameters
    unc_pars_dot = inv(R) * Y' * s;
end
