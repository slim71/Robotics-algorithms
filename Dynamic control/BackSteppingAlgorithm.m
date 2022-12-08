%% Startup

% Launches the file that prepares robot and trajectory
if ~exist('abbirb', 'var')
    control_algorithms
end

%% Simulation

% Initial conditions
% q0 defined in the preparation file
dq0 = [0, 0, 0, 0, 0, 0];

result_bs_q = zeros(size(q0, 2), length(t));
result_bs_dq = zeros(size(q0, 2), length(t));
result_bs_ddq = zeros(size(q0, 2), length(t));

Lambda = 100 * diag([1 1 1 1 1 1]);
Kd = Lambda;

% Lambda = diag([0.3 0.3 0.1 0.1 0.05 0.01]);
% Kd = 1 * diag([1 1 1 1 1 1]);
Kp = 1 * diag([1 1 1 1 1 1]);

q0s = [q0'; dq0'];
for i = 1:length(t)
    des = [q_des(:, i); dq_des(:, i); ddq_des(:, i)];
    tic
    % Output arguments: 
    % 1. time (not used)
    % 2. [q, dq], stacked in rows
    [~, bs_result] = ode15s(@(t,y) BS_odefun(t,y,des,Lambda,Kd,Kp,abbirb), ...
                            [t_init, t_end/100], q0s);
%     [~, bs_result] = ode15s(@(t,y) BS_odefun(t,y,des,Lambda,Kd,Kp,abbirb, des_pose(i,:), des_twist(i,:)), ...
%                             [t_init, t_end/100], q0s);
    fprintf("Elapsed %d s for BS solution, %d iterations \n", toc, i);

    % Save results
    result_bs_q(:, i) = bs_result(end, 1:6)';
    result_bs_dq(:, i) = bs_result(end, 7:12)';

    % Change initial conditions to the newly-reached point
    q0s = [bs_result(end, 1:6)'; bs_result(end, 7:end)'];
end

%% Plotting result data

bs_ext_forces = zeros(size(result_bs_q));
bs_acc = zeros(size(result_bs_q));
for i = 1:size(result_bs_q, 2)
    [bs_ext_forces(:, i), bs_acc(:, i)] = BS_ExtForcesAndAcc(abbirb, result_bs_q(:, i)', ...
                                                   result_bs_dq(:, i)', ...
                                                   q_des(:, i)', ...
                                                   dq_des(:, i)', ...
                                                   ddq_des(:, i)', ...
                                                   Lambda,Kd,Kp);
end

resq_fig = figure2('Name', 'Resulting BS joint configuration');
sgtitle("Resulting BS joint configuration");
n_resq = size(result_bs_q, 1);
for i = 1:n_resq
    sp = subplot(n_resq, 1, i);
    hold on
    
    plot(t, rem(result_bs_q(i, :), 2*pi))

    grid
    xlabel("time [s]");
    ylh = ylabel(qlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resdq_fig = figure2('Name', 'Resulting BS joint velocities');
sgtitle("Resulting BS joint velocities");
n_resdq = size(result_bs_dq, 1);
for i = 1:n_resdq
    sp = subplot(n_resdq, 1, i);
    hold on
    
    plot(t, result_bs_dq(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qdlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resddq_fig = figure2('Name', 'Resulting BS joint accelerations');
sgtitle("Resulting BS joint accelerations");
n_resddq = size(bs_acc, 1);
for i = 1:n_resddq
    sp = subplot(n_resddq, 1, i);
    hold on
    
    plot(t, bs_acc(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qddlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

bs_pose_error = zeros(6, size(result_bs_q, 2));
bs_pose_angles = zeros(3, size(result_bs_q, 2));
% Compute pose errors
for i = 1:size(result_bs_q, 2)
    tempfk = abbirb.fkine(result_bs_q(:, i)');
    bs_pose_error(1:3, i) = indexAt(des_pose(i, :), 1:3) - indexAt(tempfk.T, 1:3, 4)';
    bs_pose_angles(:, i) = indexAt(MatToRPY(indexAt(tempfk.T, 1:3, 1:3)), 3:-1:1);
    bs_pose_error(4:6, i) = indexAt(des_pose(i, :), 4:6) - bs_pose_angles(:, i)';
end

jposerr_fig = figure2('Name', 'BS end-effector pose errors');
sgtitle("BS end-effector pose errors");
n_poseerr = size(bs_pose_error, 1);
for i = 1:n_poseerr
    sp = subplot(n_poseerr, 1, i);
    hold on
    
    if i <= 3
        plot(t, bs_pose_error(i, :))
    else
        plot(t, wrapTo180(rad2deg(bs_pose_error(i, :))))
    end

    grid
    xlabel("time [s]");
    ylh = ylabel(eeaxislabels(i), 'Interpreter', 'latex', 'rotation', 90, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

%% Final behavior

% Plotting desired and actual EE position
bs_ee_pos = zeros(3, size(result_bs_q, 2));
for i = 1:size(result_bs_q, 2)
    fki = abbirb.fkine(result_bs_q(:, i));
    bs_ee_pos(:, i) = indexAt(fki.T, 1:3, 4);
end

bs_ee_x = bs_ee_pos(1, :);
bs_ee_y = bs_ee_pos(2, :);
bs_ee_z = bs_ee_pos(3, :);

% Be aware of the magnitude of the error, even if it seems 'a lot'
comp_fig = figure2('Name', 'Comparing result and reference');
sgtitle("Resulting BS e-e positions");
plot3(x_traj, y_traj, z_traj, '-or');
grid on
hold on
plot3(bs_ee_x, bs_ee_y, bs_ee_z, '-ob');
plot3(bs_ee_x(1), bs_ee_y(1), bs_ee_z(1), 'g*', 'MarkerSize', 20);
plot3(bs_ee_x(end), bs_ee_y(end), bs_ee_z(end), 'm*', 'MarkerSize', 20);
legend("Desired trajectory", "E-E positions", "Starting position", "Ending position");
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')

% Keeping above in mind, look at this
final_dif = figure2('Name', 'Resulting BS behavior');
for i = 1:size(result_bs_q, 2)
    abbirb.plot(result_bs_q(:, i)')
end

