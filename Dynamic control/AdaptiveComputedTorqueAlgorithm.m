%% Startup

% Launches the file that prepares robot and trajectory
if ~exist('abbirb', 'var')
    control_algorithms
end

%% Uncertainties
% Copy robot object to isolate computations
abbact = abbirb_unc;

n_links = size(abbact.links, 2);
n_link_pars = 10;

% Initial values for dynamics parameters
% Vertical stack: joints
% Horizontal stack: iterations
for i = 1:n_links
    act_unc_pars((i-1)*n_link_pars+1:i*n_link_pars, 1) = [abbact.links(i).m, ... % masses
                               abbact.links(i).r, ... % com
                               abbact.links(i).I(1,1), ... % inertia
                               abbact.links(i).I(1,2), ...
                               abbact.links(i).I(1,3), ...
                               abbact.links(i).I(2,2), ...
                               abbact.links(i).I(2,3), ...
                               abbact.links(i).I(3,3)]';
end

%% Simulation

% Initial conditions
% q0 defined in the preparation file
dq0 = [0, 0, 0, 0, 0, 0];

result_act_q = zeros(size(q0, 2), length(t));
result_act_dq = zeros(size(q0, 2), length(t));
result_act_ddq = zeros(size(q0, 2), length(t));

% Gains
kp = 100 * diag([1 1 1 1 1 1]);
kd = 100 * diag([1 1 1 1 1 1]);
kg = 100 * diag([1 1 1 1 1 1]);
% R,P in Lyapunov candidate, so they're positive definite
R = 100 * diag(ones(1, n_link_pars*n_links));
P = 100 * diag(ones(1, 2*n_links));

q0s = [q0'; dq0'];
for i = 1:length(t)
    des = [q_des(:, i); dq_des(:, i); ddq_des(:, i)];
    
    % Set dynamics parameters
    for j = 1:n_links
        % Masses
        abbact.links(j).m = act_unc_pars((j-1)*n_link_pars+1, i);
        % Centers of mass
        abbact.links(j).r = act_unc_pars((j-1)*n_link_pars+2:(j-1)*n_link_pars+4, i)';
        % Inertia matrices
        abbact.links(j).I = [
                      act_unc_pars((j-1)*n_link_pars+5:(j-1)*n_link_pars+7, i)';
                     [act_unc_pars((j-1)*n_link_pars+6, i)', ...
                      act_unc_pars((j-1)*n_link_pars+8:(j-1)*n_link_pars+9, i)'];
                     [act_unc_pars((j-1)*n_link_pars+7, i)', ...
                      act_unc_pars((j-1)*n_link_pars+9:(j-1)*n_link_pars+10, i)']
                    ];
    end
    
    
    tic
    % Output arguments: 
    % 1. time (not used)
    % 2. [q, dq], stacked in rows
    [~, act_result] = ode15s(@(t,y) CT_odefun(t,y,des,kp,kd,abbact), ...
                             [t_init, t_end], q0s);
    fprintf("Elapsed %d s for ACT solution, %d iterations \n", toc, i);

    % Save results
    result_act_q(:, i) = act_result(end, 1:6)';
    result_act_dq(:, i) = act_result(end, 7:12)';

    % Change initial conditions to the newly-reached point
    q0s = [act_result(end, 1:6)'; act_result(end, 7:end)'];

    % Update dynamics parameters
    if i < length(t)
        tic
        unc_odeset = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
        [~, int_pars] = ode15s(@(t,y) ACTUncertainties_odefun(t,y, ...
                                        [result_act_q(:,i)',result_act_dq(:,i)'], ...
                                                des,R,P,abbact), ...
                                        [t_init, t_end], act_unc_pars(:, i));
        fprintf("Elapsed %d s for ACT uncertainties, %d iterations \n", toc, i);
        act_unc_pars(:, i+1) = int_pars(end, :)';
    end
    
end

%% Plotting result data
 
% Compute external force and joint accelerations
act_ext_forces = zeros(size(result_act_q));
act_acc = zeros(size(result_act_q));
for i = 1:size(result_act_q, 2)

    % Set dynamics parameters
    for j = 1:n_links
        % Masses
        abbact.links(j).m = act_unc_pars((j-1)*n_link_pars+1, i);
        % Centers of mass
        abbact.links(j).r = act_unc_pars((j-1)*n_link_pars+2:(j-1)*n_link_pars+4, i)';
        % Inertia matrices
        abbact.links(j).I = [
                      act_unc_pars((j-1)*n_link_pars+5:(j-1)*n_link_pars+7, i)';
                     [act_unc_pars((j-1)*n_link_pars+6, i)', ...
                      act_unc_pars((j-1)*n_link_pars+8:(j-1)*n_link_pars+9, i)'];
                     [act_unc_pars((j-1)*n_link_pars+7, i)', ...
                      act_unc_pars((j-1)*n_link_pars+9:(j-1)*n_link_pars+10, i)']
                    ];
    end

    [act_ext_forces(:, i), act_acc(:, i)] = CT_ExtForcesAndAcc(abbact, result_act_q(:, i)', ...
                                                   result_act_dq(:, i)', ...
                                                   q_des(:, i)', ...
                                                   dq_des(:, i)', ...
                                                   ddq_des(:, i)', ...
                                                   kp, kd);
end

resq_fig = figure2('Name', 'Resulting ACT joint configuration');
sgtitle("Resulting ACT joint configuration");
n_resq = size(result_act_q, 1);
for i = 1:n_resq
    sp = subplot(n_resq, 1, i);
    hold on
    
    plot(t, rem(result_act_q(i, :), 2*pi))

    grid
    xlabel("time [s]");
    ylh = ylabel(qlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resdq_fig = figure2('Name', 'Resulting ACT joint velocities');
sgtitle("Resulting ACT joint velocities");
n_resdq = size(result_act_dq, 1);
for i = 1:n_resdq
    sp = subplot(n_resdq, 1, i);
    hold on
    
    plot(t, result_act_dq(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qdlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resddq_fig = figure2('Name', 'Resulting ACT joint accelerations');
sgtitle("Resulting ACT joint accelerations");
n_resddq = size(act_acc, 1);
for i = 1:n_resddq
    sp = subplot(n_resddq, 1, i);
    hold on
    
    plot(t, act_acc(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qddlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

act_pose_error = zeros(6, size(result_act_q, 2));
act_pose_angles = zeros(3, size(result_ct_q, 2));
% Compute pose errors
for i = 1:size(result_act_q, 2)
    tempfk = abbact.fkine(result_act_q(:, i)');
    act_pose_error(1:3, i) = indexAt(des_pose(i, :), 1:3) - indexAt(tempfk.T, 1:3, 4)';
    act_pose_angles(:, i) = indexAt(MatToRPY(indexAt(tempfk.T, 1:3, 1:3)), 3:-1:1);
    act_pose_error(4:6, i) = indexAt(des_pose(i, :), 4:6) - act_pose_angles(:, i)';
end

jposerr_fig = figure2('Name', 'ACT end-effector pose errors');
sgtitle("ACT end-effector pose errors");
n_poseerr = size(act_pose_error, 1);
for i = 1:n_poseerr
    sp = subplot(n_poseerr, 1, i);
    hold on
    
    if i <= 3
        plot(t, act_pose_error(i, :))
    else
        plot(t, wrapTo180(rad2deg(act_pose_error(i, :))))
    end

    grid
    xlabel("time [s]");
    ylh = ylabel(eeaxislabels(i), 'Interpreter', 'latex', 'rotation', 90, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

%% Final behavior

% Plotting desired and actual EE position
result_act_ee_pos = zeros(3, size(result_act_q, 2));
for i = 1:size(result_act_q, 2)
    fki = abbact.fkine(result_act_q(:, i));
    result_act_ee_pos(:, i) = indexAt(fki.T, 1:3, 4);
end

act_ee_x = result_act_ee_pos(1, :);
act_ee_y = result_act_ee_pos(2, :);
act_ee_z = result_act_ee_pos(3, :);

% Be aware of the magnitude of the error, even if it seems 'a lot'
comp_fig = figure2('Name', 'Comparing ACT result and reference');
sgtitle("Resulting ACT e-e positions");
plot3(x_traj, y_traj, z_traj, '-or');
grid on
hold on
plot3(act_ee_x, act_ee_y, act_ee_z, '-ob');
plot3(act_ee_x(1), act_ee_y(1), act_ee_z(1), 'g*', 'MarkerSize', 20);
plot3(act_ee_x(end), act_ee_y(end), act_ee_z(end), 'm*', 'MarkerSize', 20);
legend("Desired trajectory", "E-E positions", "Starting position", "Ending position");
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')

% Keeping above in mind, look at this
final_dif = figure2('Name', 'ACT resulting behavior');
for i = 1:size(result_act_q, 2)
    abbact.plot(result_act_q(:, i)')
end

