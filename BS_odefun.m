function ret_array = BS_odefun(t, y, des_array, gain, robot)

    n = size(y, 1)/2; % or size(des_array, 2)

    qs = y(1:n)';
    dqs = y(n+1:2*n)';

    % Simply any desired configuration derived from Cinematic Inversion
    q_des = des_array(1:n)';
    dq_des = des_array(n+1:2*n)';
    ddq_des = des_array(2*n+1:end)';
    
    M = robot.inertia(qs);
    C = robot.coriolis(qs, dqs);
    G = robot.gravload(qs);
    J = robot.jacob0(qs);

    % Errors
    e = q_des - qs;
    de = dq_des - dqs;

    dq_ref = dq_des' + gain*e';
    ddq_ref = ddq_des' + gain*de';
    s = de' + gain*e'; % dq_refs - dqs;

    % Torque
    % TODO: maybe only '+err' is enough? maybe even something else?
    tau = M*ddq_ref + C*dq_ref + G' +  gain*s + e';

    % z = dq/dt
    z = dqs;
    % dz/dt = ddq/ddt
    dzdt = pinv(M) * (tau - C*dqs' - G');
    
    % qdot, qddot
    ret_array = [z'; dzdt];
end
