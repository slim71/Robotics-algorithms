%% Startup

% Launches the file that prepares robot and trajectory
if ~exist('abbirb', 'var')
    control_algorithms
end

%% Uncertainties
% Copy robot object to isolate computations
abbabs = abbirb;

n_links = size(abbabs.links, 2);
n_link_pars = 10;

% Initial values for dynamics parameters
% Vertical stack: joints
% Horizontal stack: iterations
for i = 1:n_links
    abs_unc_pars((i-1)*n_link_pars+1:i*n_link_pars, 1) = [abbabs.links(i).m, ... % masses
                               abbabs.links(i).r, ... % com
                               abbabs.links(i).I(1,1), ... % inertia
                               abbabs.links(i).I(1,2), ...
                               abbabs.links(i).I(1,3), ...
                               abbabs.links(i).I(2,2), ...
                               abbabs.links(i).I(2,3), ...
                               abbabs.links(i).I(3,3)]';
end

%% Simulation

% Initial conditions
% q0 defined in the preparation file
dq0 = [0, 0, 0, 0, 0, 0];

result_abs_q = zeros(size(q0, 2), length(t));
result_abs_dq = zeros(size(q0, 2), length(t));
result_bs_ddq = zeros(size(q0, 2), length(t));

% Gains
Kp = 1000 * diag([1 1 1 1 1 1]);
Kd = 1000 * diag([1 1 1 1 1 1]);
Kg = 1000 * diag([1 1 1 1 1 1]);
% R in Lyapunov candidate, so it's positive definite
R = 10 * diag(ones(1, n_link_pars*n_links));

q0s = [q0'; dq0'];
for i = 1:length(t)
    des = [q_des(:, i); dq_des(:, i); ddq_des(:, i)];
    
    % Set dynamics parameters
    for j = 1:n_links
        % Masses
        abbabs.links(j).m = abs_unc_pars((j-1)*n_link_pars+1, i);
        % Centers of mass
        abbabs.links(j).r = abs_unc_pars((j-1)*n_link_pars+2:(j-1)*n_link_pars+4, i)';
        % Inertia matrices
        abbabs.links(j).I = [
                      abs_unc_pars((j-1)*n_link_pars+5:(j-1)*n_link_pars+7, i)';
                     [abs_unc_pars((j-1)*n_link_pars+6, i)', ...
                      abs_unc_pars((j-1)*n_link_pars+8:(j-1)*n_link_pars+9, i)'];
                     [abs_unc_pars((j-1)*n_link_pars+7, i)', ...
                      abs_unc_pars((j-1)*n_link_pars+9:(j-1)*n_link_pars+10, i)']
                    ];
    end
    
    tic
    % Output arguments: 
    % 1. time (not used)
    % 2. [q, dq], stacked in rows
    [~, abs_result] = ode15s(@(t,y) BS_odefun(t,y,des,Kg,Kd,Kp,abbabs), ...
                             [t_init, t_end/100], q0s);
    fprintf("Elapsed %d s for ABS solution, %d iterations \n", toc, i);

    % Save results
    result_abs_q(:, i) = abs_result(end, 1:6)';
    result_abs_dq(:, i) = abs_result(end, 7:12)';

    % Change initial conditions to the newly-reached point
    q0s = [abs_result(end, 1:6)'; abs_result(end, 7:end)'];

    % Update dynamics parameters
    if i < length(t)+1
        tic
        unc_odeset = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
        [~, int_pars] = ode15s(@(t,y) ABSUncertainties_odefun(t,y, ...
                                        [result_abs_q(:,i)',result_abs_dq(:,i)'], ...
                                        des,Kg,R), ...
                                        [t_init, t_end/100], abs_unc_pars(:, i));
        fprintf("Elapsed %d s for ABS uncertainties, %d iterations \n", toc, i);
        abs_unc_pars(:, i+1) = int_pars(end, :)';
    end
    
end

%% Plotting result data

% Compute external force and joint accelerations
abs_ext_forces = zeros(size(result_abs_q));
abs_acc = zeros(size(result_abs_q));
for i = 1:size(result_abs_q, 2)

    % Set dynamics parameters
    for j = 1:n_links
        % Masses
        abbabs.links(j).m = abs_unc_pars((j-1)*n_link_pars+1, i);
        % Centers of mass
        abbabs.links(j).r = abs_unc_pars((j-1)*n_link_pars+2:(j-1)*n_link_pars+4, i)';
        % Inertia matrices
        abbabs.links(j).I = [
                      abs_unc_pars((j-1)*n_link_pars+5:(j-1)*n_link_pars+7, i)';
                     [abs_unc_pars((j-1)*n_link_pars+6, i)', ...
                      abs_unc_pars((j-1)*n_link_pars+8:(j-1)*n_link_pars+9, i)'];
                     [abs_unc_pars((j-1)*n_link_pars+7, i)', ...
                      abs_unc_pars((j-1)*n_link_pars+9:(j-1)*n_link_pars+10, i)']
                    ];
    end

    % Compute forces and accelerations
    [abs_ext_forces(:, i), abs_acc(:, i)] = BS_ExtForcesAndAcc(abbabs, ...
                                                   result_abs_q(:, i)', ...
                                                   result_abs_dq(:, i)', ...
                                                   q_des(:, i)', ...
                                                   dq_des(:, i)', ...
                                                   ddq_des(:, i)', ...
                                                   Kg, Kd, Kp);
end

resq_fig = figure2('Name', 'Resulting ABS joint configuration');
sgtitle("Resulting ABS joint configuration");
n_resq = size(result_abs_q, 1);
for i = 1:n_resq
    sp = subplot(n_resq, 1, i);
    hold on
    
    plot(t, result_abs_q(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resdq_fig = figure2('Name', 'Resulting ABS joint velocities');
sgtitle("Resulting ABS joint velocities");
n_resdq = size(result_abs_dq, 1);
for i = 1:n_resdq
    sp = subplot(n_resdq, 1, i);
    hold on
    
    plot(t, result_abs_dq(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qdlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

resddq_fig = figure2('Name', 'Resulting ABS joint accelerations');
sgtitle("Resulting ABS joint accelerations");
n_resddq = size(abs_acc, 1);
for i = 1:n_resddq
    sp = subplot(n_resddq, 1, i);
    hold on
    
    plot(t, abs_acc(i, :))

    grid
    xlabel("time [s]");
    ylh = ylabel(qddlabels(i), 'Interpreter', 'latex', 'rotation', 0, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

abs_pose_error = zeros(6, size(result_abs_q, 2));
% Compute pose errors
for i = 1:size(result_abs_q, 2)
    tempfk = abbabs.fkine(result_abs_q(:, i)');
    abs_pose_error(1:3, i) = indexAt(tempfk.T, 1:3, 4)' - [x_traj(i), y_traj(i), z_traj(i)];
    abs_pose_error(4:6, i) = MatToRPY(indexAt(tempfk.T, 1:3, 1:3) - [des_theta(i), des_phi(i), des_psi(i)]);
end

jposerr_fig = figure2('Name', 'ABS end-effector pose errors');
sgtitle("ABS end-effector pose errors");
n_poseerr = size(abs_pose_error, 1);
for i = 1:n_poseerr
    sp = subplot(n_poseerr, 1, i);
    hold on
    
    if i <= 3
        plot(t, abs_pose_error(i, :))
    else
        plot(t, rad2deg(abs_pose_error(i, :)))
    end

    grid
    xlabel("time [s]");
    ylh = ylabel(eeaxislabels(i), 'Interpreter', 'latex', 'rotation', 90, 'VerticalAlignment', 'middle');
    % Half the y-axis, left of the x-axis
    ylh.Position = [sp.XLim(1)-0.3, sum(sp.YLim)/2];
end

%% Final behavior

% Plotting desired and actual EE position
abs_ee_pos = zeros(3, size(result_abs_q, 2));
for i = 1:size(result_abs_q, 2)
    fki = abbabs.fkine(result_abs_q(:, i));
    abs_ee_pos(:, i) = indexAt(fki.T, 1:3, 4);
end

abs_ee_x = abs_ee_pos(1, :);
abs_ee_y = abs_ee_pos(2, :);
abs_ee_z = abs_ee_pos(3, :);

% Be aware of the magnitude of the error, even if it seems 'a lot'
comp_fig = figure2('Name', 'Comparing ABS result and reference');
sgtitle("Resulting ABS e-e positions");
plot3(x_traj, y_traj, z_traj, '-or');
grid on
hold on
plot3(abs_ee_x, abs_ee_y, abs_ee_z, '-ob');
xlabel('x')
ylabel('y')
zlabel('z')

% Keeping above in mind, look at this
final_dif = figure2('Name', 'ABS resulting behavior');
for i = 1:size(result_abs_q, 2)
    abbabs.plot(result_abs_q(:, i)')
end
