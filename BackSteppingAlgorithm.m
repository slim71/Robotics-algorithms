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

Kd = 10 * diag([1 1 1 1 1 1]);

q0s = [q0'; dq0'];
for i = 1:length(t)+1
    des = [q_des(:, i); dq_des(:, i); ddq_des(:, i)];
    tic
    % Output arguments: 
    % 1. time (not used)
    % 2. [q, dq], stacked in rows
    [~, bs_result] = ode15s(@(t,y) BS_odefun(t,y,des,Kd,abbirb), [t_init, t_end/100], q0s);
    fprintf("Elapsed %d s, %d iterations \n", toc, i);

    % Save results
    result_bs_q(:, i) = bs_result(end, 1:6)';
    result_bs_dq(:, i) = bs_result(end, 7:12)';

    % Change initial conditions to the newly-reached point
    q0s = [bs_result(end, 1:6)'; bs_result(end, 7:end)'];
end

%% Plotting results data

resq_fig = figure2('Name', 'Resulting joint configuration');
sgtitle("Resulting joint configuration");
n_resq = size(result_bs_q, 1);
for i = 1:n_resq
    resq_sp = subplot(n_resq, 1, i);
    hold on
    
    plot(t, result_bs_q(i, 2:end))

    grid
    xlabel("time [s]");
    ylh = ylabel(qlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    ylh.Position = [-0.3, 0, -1];
end

resdq_fig = figure2('Name', 'Resulting joint velocities');
sgtitle("Resulting joint velocities");
n_resdq = size(result_bs_dq, 1);
for i = 1:n_resdq
    resdq_sp = subplot(n_resdq, 1, i);
    hold on
    
    plot(t, result_bs_dq(i, 2:end))

    grid
    xlabel("time [s]");
    ylh = ylabel(qdlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    ylh.Position = [-0.3, 0, -1];
end

resddq_fig = figure2('Name', 'Resulting joint accelerations');
sgtitle("Resulting joint accelerations");
n_resddq = size(result_bs_ddq, 1);
for i = 1:n_resddq
    resq_sp = subplot(n_resddq, 1, i);
    hold on
    
    plot(t, result_bs_ddq(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qddlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    ylh.Position = [-0.3, 0, -1];
end

ext_forces = zeros(size(result_bs_q));
acc = zeros(size(result_bs_q));
for i = 1:size(result_bs_q, 2)
    [ext_forces(:, i), acc(:, i)] = BackStepping(abbirb, result_bs_q(:, i)', ...
                                                   result_bs_dq(:, i)', ...
                                                   q_des(:, i)', ...
                                                   dq_des(:, i)', ...
                                                   ddq_des(:, i)', ...
                                                   Kd);
end

bs_ee_x = bs_ee_pos(1, :);
bs_ee_y = bs_ee_pos(2, :);
bs_ee_z = bs_ee_pos(3, :);

% Plotting desired and actual EE position
bs_ee_pos = zeros(3, size(result_bs_q, 2));
for i = 1:size(result_bs_q, 2)
    fki = abbirb.fkine(result_bs_q(:, i));
    bs_ee_pos(:, i) = indexAt(fki.T, 1:3, 4);
end

% Be aware of the magnitude of the error, even if it seems 'a lot'
comp_fig = figure2('Name', 'Comparing result and reference');
plot3(x_traj, y_traj, z_traj, '-or');
grid on
hold on
plot3(bs_ee_x, bs_ee_y, bs_ee_z, '-ob');
xlabel('x')
ylabel('y')
zlabel('z')

% Keeping above in mind, look at this
final_dif = figure2('Name', 'Resulting behavior');
for i = 1:size(result_bs_q, 2)
    abbirb.plot(result_bs_q(:, i)')
end

