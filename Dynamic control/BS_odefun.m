function ret_array = BS_odefun(t, y, des_array, lambda, kd, kp, robot, pose_des, twist_des)

    n = size(y, 1)/2; % or size(des_array, 2)

    % Current joint configuration
    qs = y(1:n)';
    dqs = y(n+1:2*n)';

    % Desired joint configuration
    qds = des_array(1:n)';
    dqds = des_array(n+1:2*n)';
    ddqds = des_array(2*n+1:end)';

    % Robot dynamics matrices
    M = robot.inertia(qs); % Inertia
    C = robot.coriolis(qs, dqs); % Coriolis
    G = robot.gravload(qs); % Gravity
%     J = robot.jacob0(qs);

    % Errors
    e = (qds - qs);
    de = (dqds - dqs);
%     current_hom = robot.fkine(qs);
%     current_pose = [indexAt(current_hom.T, 1:3, 4)', indexAt(MatToRPY(indexAt(current_hom.T,1:3,1:3)), 3:-1:1)];
%     desired_pose = pose_des;
%     e = desired_pose - current_pose;
%     de = diff([e; zeros(1,6)],1,2)/0.1;

    dq_ref = dqds' - lambda*e';
    ddq_ref = ddqds' - lambda*de'; %derivative of dq_ref
%     dq_ref = twist_des' + lambda * e';
%     ddq_ref = diff([dq_ref, zeros(6,1)],1,2)/0.1;

    s = dq_ref - dqs';

    % Control torque
    tau_barred = e';
%     tau_barred = J * e';
    tau = M*ddq_ref + C*dq_ref + G' + kd*s + tau_barred;

    % ODE outputs
    % z = dq/dt
    z = dqs;
    % dz/dt = ddq/ddt
    dzdt = pinv(M) * (tau - C*dqs' - G');
    
    % [qdot; qddot]
    ret_array = [z'; dzdt];
end
