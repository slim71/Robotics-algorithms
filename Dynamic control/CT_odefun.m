function ret_array = CT_odefun(t, y, des_array, Kp, Kd, robot)

    n = size(y, 1)/2; % or size(des_array, 2)

    qs = y(1:n)';
    dqs = y(n+1:2*n)';

    qds = des_array(1:n)';
    dqds = des_array(n+1:2*n)';
    ddqds = des_array(2*n+1:end);
    
    % Robot dynamics matrices
    M = robot.inertia(qs);
    C = robot.coriolis(qs, dqs);
    G = robot.gravload(qs);

    % Errors
    e = (qds - qs);
    de = (dqds - dqs);

    tau = M*(ddqds + Kd*de' + Kp*e') + C*dqs' + G';

    % z = dq/dt
    z = dqs;
    % dz/dt = ddq/ddt
    dzdt = pinv(M) * (tau - C*dqs' - G');
    
    % qdot, qddot
    ret_array = [z'; dzdt];
end