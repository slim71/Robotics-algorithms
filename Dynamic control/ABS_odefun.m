function ret_array = ABS_odefun(t, y, des_array, lambda, kd, kp, robot)
% TODO: delete if BS only is ok

    n = size(y, 1)/2; % or size(des_array, 2)

    qs = y(1:n)';
    dqs = y(n+1:2*n)';

    % Simply any desired configuration derived from Cinematic Inversion
    q_des = des_array(1:n)';
    dq_des = des_array(n+1:2*n)';
    ddq_des = des_array(2*n+1:end)';
    
    % Errors
    e = q_des - qs;
    de = dq_des - dqs;

    dq_ref = dq_des' + lambda*e';
    ddq_ref = ddq_des' + lambda*de';
    s = de' + lambda*e';

    % Compute dynamics matrices
    M = robot.inertia(qs); % inertia
    C = robot.coriolis(qs, dqs); % Coriolis
    G = robot.gravload(qs); % Gravity

    tau = M*ddq_ref + C*dq_ref + G' + kd*s + kp*e';

    % z = dq/dt
    z = dqs;
    % dz/dt = ddq/ddt
    dzdt = pinv(M) * (tau - C*dqs' - G');
    
    % qdot, qddot
    ret_array = [z'; dzdt];
end
