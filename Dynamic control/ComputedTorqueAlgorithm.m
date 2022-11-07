%% Startup

% Launches the file that prepares robot and trajectory
if ~exist('abbirb', 'var')
    control_algorithms
end

%% Simulation

% Initial conditions
% q0 defined in the preparation file
dq0 = [0, 0, 0, 0, 0, 0];

result_ct_q = zeros(size(q0, 2), length(t));
result_ct_dq = zeros(size(q0, 2), length(t));
result_ct_ddq = zeros(size(q0, 2), length(t));

kp = 20 * diag([1 1 1 1 1 1]);
kd = 10 * diag([1 1 1 1 1 1]);

q0s = [q0'; dq0'];
for i = 1:length(t)
    des = [q_des(:, i); dq_des(:, i); ddq_des(:, i)];
    tic
    % Output arguments: 
    % 1. time (not used)
    % 2. [q, dq], stacked in rows
    [~, comp_result] = ode15s(@(t,y) CT_odefun(t,y,des,kp,kd,abbirb), ...
                              [t_init, t_end/100], q0s);
    fprintf("Elapsed %d s for CT solution, %d iterations \n", toc, i);

    % Save results
    result_ct_q(:, i) = comp_result(end, 1:6)';
    result_ct_dq(:, i) = comp_result(end, 7:12)';

    % Change initial conditions to the newly-reached point
    q0s = [comp_result(end, 1:6)'; comp_result(end, 7:end)'];
end

%% Plotting result data

ct_ext_forces = zeros(size(result_ct_q));
ct_acc = zeros(size(result_ct_q));
for i = 1:size(result_ct_q, 2)
    [ct_ext_forces(:, i), ct_acc(:, i)] = CT_ExtForcesAndAcc(abbirb, result_ct_q(:, i)', ...
                                                   result_ct_dq(:, i)', ...
                                                   q_des(:, i)', ...
                                                   dq_des(:, i)', ...
                                                   ddq_des(:, i)', ...
                                                   kp, kd);
end

resq_fig = figure2('Name', 'Resulting CT joint configuration');
sgtitle("Resulting CT joint configuration");
n_resq = size(result_ct_q, 1);
for i = 1:n_resq
    sp = subplot(n_resq, 1, i);
    hold on
    
    plot(t, rem(result_ct_q(i, :), 2*pi))

    grid
    xlabel("time [s]");
    ylh = ylabel(qlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resdq_fig = figure2('Name', 'Resulting CT joint velocities');
sgtitle("Resulting CT joint velocities");
n_resdq = size(result_ct_dq, 1);
for i = 1:n_resdq
    sp = subplot(n_resdq, 1, i);
    hold on
    
    plot(t, result_ct_dq(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qdlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resddq_fig = figure2('Name', 'Resulting CT joint accelerations');
sgtitle("Resulting CT joint accelerations");
n_resddq = size(ct_acc, 1);
for i = 1:n_resddq
    sp = subplot(n_resddq, 1, i);
    hold on
    
    plot(t, ct_acc(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qddlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

ct_pose_error = zeros(6, size(result_ct_q, 2));
% Compute pose errors
for i = 1:size(result_ct_q, 2)
    tempfk = abbirb.fkine(result_ct_q(:, i)');
    ct_pose_error(1:3, i) = indexAt(tempfk.T, 1:3, 4)' - [x_traj(i), y_traj(i), z_traj(i)];
    ct_pose_error(4:6, i) = MatToRPY(indexAt(tempfk.T, 1:3, 1:3) - [des_theta(i), des_phi(i), des_psi(i)]);
end

jposerr_fig = figure2('Name', 'CT end-effector pose errors');
sgtitle("CT end-effector pose errors");
n_poseerr = size(ct_pose_error, 1);
for i = 1:n_poseerr
    sp = subplot(n_poseerr, 1, i);
    hold on
    
    if i <= 3
        plot(t, ct_pose_error(i, :))
    else
        plot(t, wrapTo180(rad2deg(ct_pose_error(i, :))))
    end

    grid
    xlabel("time [s]");
    ylh = ylabel(eeaxislabels(i), 'Interpreter', 'latex', 'rotation', 90, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

%% Final behavior

% Plotting desired and actual EE position
result_ee_pos = zeros(3, size(result_ct_q, 2));
for i = 1:size(result_ct_q, 2)
    fki = abbirb.fkine(result_ct_q(:, i));
    result_ee_pos(:, i) = indexAt(fki.T, 1:3, 4);
end

ee_x = result_ee_pos(1, :);
ee_y = result_ee_pos(2, :);
ee_z = result_ee_pos(3, :);

% Be aware of the magnitude of the error, even if it seems 'a lot'
comp_fig = figure2('Name', 'Comparing result and reference');
sgtitle("Resulting CT e-e positions");
plot3(x_traj, y_traj, z_traj, '-or');
grid on
hold on
plot3(ee_x, ee_y, ee_z, '-ob');
plot3(ee_x(1), ee_y(1), ee_z(1), 'g*', 'MarkerSize', 20);
plot3(ee_x(end), ee_y(end), ee_z(end), 'm*', 'MarkerSize', 20);
legend("Desired trajectory", "E-E positions", "Starting position", "Ending position");
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')

% Keeping above in mind, look at this
final_dif = figure2('Name', 'Resulting CT behavior');
for i = 1:size(result_ct_q, 2)
    abbirb.plot(result_ct_q(:, i)')
end

