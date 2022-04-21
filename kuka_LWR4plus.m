%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data [cm]
l0 = 0.11; l1 = 0.20; l2 = 0.20; l3 = 0.20; l4 = 0.20; l5 = 0.19;

%% Denavit-Hartemberg and Jacobian
%                a   alpha   d       theta    joint_type
denavit = @(q)([[0,  pi/2,   0,      q(1),       "R"] 
                [0,  -pi/2,  0,      q(2),       "R"] 
                [0,  -pi/2,  l2+l3,  q(3),       "R"] 
                [0,  pi/2,   0,      q(4),       "R"] 
                [0,  pi/2,   l4+l5,  q(5),       "R"] 
                [0,  -pi/2,  0,      q(6),       "R"] 
                [0,  0,      0,      q(7),       "R"]]);

J = @(q) DHJacob0(denavit(q), 0);
J4 = @(q) DHJacob0(denavit(q), 4);

ForKine = @(q) DHFKine(denavit(q));
FK4 = @(q) DHFKine(denavit(q), 4);

%% Tasks
t_task = 0:pi/10:3*pi;

% Desired EE motion
despose = @(t) [0.1 * sin(t), 0.2 * cos(t), 0.65, t, t, t];
% Make actual derivative of function, avoiding 'diff' problems
destwist = matlabFunction(diff(despose(t)));

% As function of t and q only for consistency
despose4 = @(t) [0 0 0 0 0 0];
destwist4 = @(t,q) [0 0 0 0 0 0];

traj_fig = figure2('Name', 'Task trajectory');
title("Trajectory");
grid
hold on
first_traj = [];
for k = t_task
    first_traj = [first_traj; despose(k)];
end
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-')
view(45, 45)

%% Robot design
% TODO: build upon denavit
kuka = rigidBodyTree('Dataformat', 'column');

% It's sufficient to drop the last column and give whatever values to q:
% joint variables are ignored in setFixedTransform
dhparams = double(denavit(zeros(7, 1)));
dhparams = dhparams(:, 1:end-1);

bodyNames = {'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'};
parentNames = {'base', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6'};
jointNames = {'j1', 'j2', 'j3', 'j4', 'j5', 'j6', 'j7'};
jointTypes = {'revolute', 'revolute', 'revolute', 'revolute', 'revolute', 'revolute', 'revolute'};

for k = 1:size(bodyNames, 2)
    % Create a rigidBody object with a unique name
    kukaBodies(k) = rigidBody(bodyNames{k});

    % Create a rigidBodyJoint object and give it a unique name
    kukaBodies(k).Joint = rigidBodyJoint(jointNames{k}, jointTypes{k});

    % Use setFixedTransform to specify the body-to-body transformation using DH parameters
    setFixedTransform(kukaBodies(k).Joint, dhparams(k,:), 'dh');

    % Attach the body joint to the robot
    addBody(kuka, kukaBodies(k), parentNames{k});
end

showdetails(kuka)

%% Needed functions
curr_pos = @(q) indexAt(ForKine(q), 1:3, 4);
des_pos = @(t) indexAt(despose(t), 1:3);
pose_error = @(t, q) (des_pos(t) - curr_pos(q)');
J4_pos = @(q) indexAt(FK4(q), 1:3, 4);
J4_despos = @(t) indexAt(despose4(t), 1:3);
J4_pose_error = @(t, q) (J4_despos(t) - J4_pos(q)');

curr_or = @(q) MatToQuat(indexAt(ForKine(q),1:3,1:3));
des_or = @(t) EulToQuat(indexAt(despose(t),4:6), 'ZYX', false);
or_error = @(t, q) indexAt(QuatProd(des_or(t), QuatInv(curr_or(q))), 2:4);
j4_or = @(q) MatToQuat(indexAt(FK4(q),1:3,1:3));
J4_desor = @(t) EulToQuat(indexAt(despose4(t),4:6), 'ZYX', false);
J4_orerr = @(t, q) indexAt(QuatProd(J4_desor(t), QuatInv(j4_or(q))), 2:4);

tot_errors = @(t, q) [pose_error(t,q), or_error(t,q)];
J4_tot_errors = @(t, q) [J4_pose_error(t,q), J4_orerr(t,q)];

%% IK parameters
qinit = [-pi/2, pi/4, 0, pi/3,0,0,0];

t0 = 0;
tf = 7;

%% SOT (1 task)
des_pos = @(t) despose(t);
curr_pos = @(q) [indexAt(ForKine(q), 1:3, 4)', ...
                    MatToRPY(indexAt(ForKine(q), 1:3, 1:3))];

pose_error = @(t, q) curr_pos(q) - des_pos(t);

% (-50; 0), but -10 is better
K = -10 * eye(6); %-10 and 15s seems good-ish, lambda -1
task1 = @(t, q) (K * pose_error(t, q)')';

% -3 -3 "works" but it's not accurate
% -2 -3 looks good [5s]
% -2 -4 too
% -2 -5, -3 -4 stops for tolerance
options = odeset('RelTol', 1e-2, 'AbsTol', 1e-4);

% better with -1; -5 is horrible
[t1_sot, y1_sot] = ode15s(@(t, y) sot(t, y, {J}, {task1}, 10^(-1)),...
                       [t0 tf], qinit, options);

q1_sot = y1_sot';

sot1_fig = figure2('Name', 'SoT (one task) with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka, qinit');
grid
for k = 1:size(q1_sot, 2)
    fk = ForKine(q1_sot(:, k));
    sot1_positions = fk(1:3, 4);
    plot3(sot1_positions(1), sot1_positions(2), sot1_positions(3), 'o');
end

% Animated plot
% sot1_fig2 = figure2('Name', 'SoT (one task) with ode15s');
% title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
% for i=[1:size(q1_sot, 2)]
%     plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black')
%     hold on
%     grid
%     show(kuka, q1_sot(:, i));
%     view(45, 45)
%     drawnow
%     pause(0.1)
%     hold off
% end
%% SoT (1 task) errors

sot1_err = zeros(6, size(q1_sot, 2));
for k = 1:size(q1_sot, 2)
    sot1_err(:, k) = tot_errors(t1_sot(k), q1_sot(:, k));
end

sot1err_fig = figure2('Name', 'SoT (1 task) errors');
for i = 1:6
    sp_sot1 = subplot(6, 1, i);

    if i <= 3 % position errors
        plot(t1_sot, sot1_err(i, :))
        ylabel("error [m]");

    else % orientation error
        plot(t1_sot, rad2deg(rem(sot1_err(i, :), 2*pi)))

        ylabel("error [deg]");
    end

    grid
    xlabel("time [s]");
    switch i
        case 1
            title("X-axis position");
        case 2
            title("Y-axis position");
        case 3
            title("Z-axis position");
        case 4
            title("X-axis orientation");
        case 5
            title("Y-axis orientation");
        case 6
            title("Z-axis orientation");
    end
end

%% SoT
[t_sot, y_sot] = ode15s(@(t,y) sot(t, y, {J, J4}, {task1, destwist4}, 10^(-1)),...
                       [t0 tf], qinit, options);

q_sot = y_sot';

sot_fig = figure2('Name', 'SoT with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka, qinit');
grid
for k = 1:size(q_sot, 2)
    fk = ForKine(q_sot(:, k));
    sot_positions = fk(1:3, 4);
    plot3(sot_positions(1), sot_positions(2), sot_positions(3), 'o');
end

% Animated plot
% sot2_fig = figure2('Name', 'SoT with ode15s');
% title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
% for i=[1:size(q_sot, 2)]
%     plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black')
%     hold on
%     grid
%     show(kuka, q_sot(:, i));
%     view(45, 45)
%     drawnow
%     pause(0.1)
%     hold off
% end

%% SoT errors

sot_err = zeros(6, size(q_sot, 2));
sot_err4 = zeros(6, size(q_sot, 2));
for k = 1:size(q_sot, 2)
    sot_err(:, k) = tot_errors(t_sot(k), q_sot(:, k));
    sot_err4(:, k) = J4_tot_errors(t_sot(k), q_sot(:, k));
end

soterr_fig = figure2('Name', 'SoT errors');
for i = 1:6
    sp_sot = subplot(6, 1, i);
    hold on

    if i <= 3 % position errors
        plot(t_sot, sot_err(i, :))
        plot(t_sot, sot_err4(i, :))
        ylabel("error [m]");

    else % orientation error
        plot(t_sot, rad2deg(rem(sot_err(i, :), 2*pi)))
        plot(t_sot, rad2deg(rem(sot_err4(i, :), 2*pi)))

        ylabel("error [deg]");
    end

    grid
    legend("EE", "J4");
    xlabel("time [s]");
    switch i
        case 1
            title("X-axis position");
        case 2
            title("Y-axis position");
        case 3
            title("Z-axis position");
        case 4
            title("X-axis orientation");
        case 5
            title("Y-axis orientation");
        case 6
            title("Z-axis orientation");
    end
end

%% Reverse Priority 1 task (does not complete)
% experimented with different combinations withouth changing the desired
% task first. In pratical times, these do not go on computing:
% -6 -6, -6 -5, -6 -4, -6 -3, -6 -2, -6 -1
% -5 -6, -5 -5, -5 -4, -5 -3, -5 -2, -5 -1
% -4 -6, -4 -5, -4 -4, -4 -3, -4 -2, -4 -1
% -3 -6, -3 -5, -3 -4, -3 -3, -3 -2, -3 -1
% -2 -6, -2 -5, -2 -4, -2 -3, -2 -2, -2 -1
% -1 -6, -1 -5, -1 -4, -1 -3, -1 -2, -1 -1  <- horrible
% options2 = odeset('RelTol', 1e-5, 'AbsTol', 1e-1);
%{
% starting with 10 as above
% tried but too slow in practice, so stopped:
% -1, -100, -1000, 1, 10, 100, 1000
K = -10 * eye(6); %
task2 = @(t, q) (K * pose_error(t, q)')';

[t1_rp, y1_rp] = ode45(@(t,y) rp(t, y, {J, J4}, {task2}, 10^(-2)),...
                       [t0 tf], qinit, options);

q1_rp = y1_rp';

rp1_fig = figure2('Name', 'RP with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka, qinit');
grid
for k=[1:size(q1_rp, 2)]
    fk = ForKine(q1_rp(:, k));
    rp1_position = fk(1:3, 4);
    plot3(rp1_position(1), rp1_position(2), rp1_position(3), 'o');
end

% Animated plot
% rp1_fig2 = figure2('Name', 'RP with ode15s');
% title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
% for i=[1:size(q1_rp, 2)]
%     plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black')
%     hold on
%     grid
%     show(kuka, q1_rp(:, i));
%     view(45, 45)
%     drawnow
%     pause(0.1)
%     hold off
% end
%}
%% RP (1 task) errors
% N/A since this integration does not complete

%% Reverse Priority

% with 1e-1 the errors have spikes
[t_rp, y_rp] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, destwist4}, 10^(-2)),...
                       [t0 tf], qinit, options);
[t_rpAlt, y_rpAlt] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, destwist4}, 10^(-1)),...
                       [t0 tf], qinit, options);

q_rp = y_rp';
q_rpAlt = y_rpAlt';

rp_fig = figure2('Name', 'RP with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka, qinit');
grid
for k=[1:size(q_rp, 2)]
    fk = ForKine(q_rp(:, k));
    rp_position = fk(1:3, 4);
    plot3(rp_position(1), rp_position(2), rp_position(3), 'o');
end

% Animated plot
% rp_fig2 = figure2('Name', 'RP with ode15s');
% title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
% for i=[1:size(q_rp, 2)]
%     plot3(first_traj(:, 1), first_traj(:, 2),first_traj(:, 3), '-', 'Color', 'black')
%     hold on
%     grid
%     show(kuka, q_rp(:, i));
%     view(45, 45)
%     drawnow
%     pause(0.1)
%     hold off
% end

%% RP errors

rp_err = zeros(6, size(q_rp, 2));
rp_err4 = zeros(6, size(q_rp, 2));
for k = 1:size(q_rp, 2)
    rp_err(:, k) = tot_errors(t_rp(k), q_rp(:, k));
    rp_err4(:, k) = J4_tot_errors(t_rp(k), q_rp(:, k));
end

rp_errAlt = zeros(6, size(q_rpAlt, 2));
rp_errAlt4 = zeros(6, size(q_rpAlt, 2));
for k = 1:size(q_rpAlt, 2)
    rp_errAlt(:, k) = tot_errors(t_rpAlt(k), q_rpAlt(:, k));
    rp_errAlt4(:, k) = J4_tot_errors(t_rpAlt(k), q_rpAlt(:, k));
end

rperr_fig = figure2('Name', 'RP errors');
for i = 1:6
    sp_rp = subplot(6, 1, i);
    hold on

    if i <= 3 % position errors
        plot(t_rp, rp_err(i, :))
        plot(t_rp, rp_err4(i, :))

        ylabel("error [m]");

    else % orientation error
        plot(t_rp, rad2deg(rem(rp_err(i, :), 2*pi)))
        plot(t_rp, rad2deg(rem(rp_err4(i, :), 2*pi)))

        ylabel("error [deg]");

    end

    grid
    legend("EE", "J4");
    xlabel("time [s]");
    switch i
        case 1
            title("X-axis position");
        case 2
            title("Y-axis position");
        case 3
            title("Z-axis position");
        case 4
            title("X-axis orientation");
        case 5
            title("Y-axis orientation");
        case 6
            title("Z-axis orientation");
    end
end

rpAlterr_fig = figure2('Name', 'RP (alternate version) errors');
for i = 1:6
    sp_rp = subplot(6, 1, i);
    hold on

    if i <= 3 % position errors
        plot(t_rpAlt, rp_errAlt(i, :))
        plot(t_rpAlt, rp_errAlt4(i, :))

        ylabel("error [m]");

    else % orientation error
        plot(t_rpAlt, rad2deg(rem(rp_errAlt(i, :), 2*pi)))
        plot(t_rpAlt, rad2deg(rem(rp_errAlt4(i, :), 2*pi)))

        ylabel("error [deg]");

    end

    grid
    legend("EE", "J4");
    xlabel("time [s]");
    switch i
        case 1
            title("X-axis position");
        case 2
            title("Y-axis position");
        case 3
            title("Z-axis position");
        case 4
            title("X-axis orientation");
        case 5
            title("Y-axis orientation");
        case 6
            title("Z-axis orientation");
    end
end

%% E-E error comparison

poserr_fig = figure2('Name', 'EE position errors comparison');
for i = 1:3
    subplot(3, 1, i)
    hold on

    plot(t1_sot, sot1_err(i, :))
    plot(t_sot, sot_err(i, :))
    plot(t_rp, rp_err(i, :))
    plot(t_rpAlt, rp_errAlt(i, :))

    grid
    legend("SoT (1 task)", "SoT", "RP", "RP (1e-1)");
    xlabel("time [s]");
    ylabel("error [m]");
    switch i
        case 1
            title("X-axis");
        case 2
            title("Y-axis");
        case 3
            title("Z-axis");
    end
end

orerr_fig = figure2('Name', 'EE orientation errors comparison');
for i = 4:6
    sp_orerr = subplot(3, 1, i-3);
    hold on

    plot(t1_sot, rad2deg(rem(sot1_err(i, :), 2*pi)))
    plot(t_sot, rad2deg(rem(sot_err(i, :), 2*pi)))
    plot(t_rp, rad2deg(rem(rp_err(i, :), 2*pi)))
    plot(t_rpAlt, rad2deg(rem(rp_errAlt(i, :), 2*pi)))

    grid
    legend("SoT (1 task)", "SoT", "RP", "RP (1e-1)");
    xlabel("time [s]");
    ylabel("error [deg]");
    switch i
        case 1
            title("X-axis");
        case 2
            title("Y-axis");
        case 3
            title("Z-axis");
    end
end

%% 4th joint error comparison

J4_poserr_fig = figure2('Name', '4th joint position errors comparison');
for i = 1:3
    subplot(3, 1, i)
    hold on

    plot(t_sot, sot_err4(i, :))
    plot(t_rp, rp_err4(i, :))
    plot(t_rpAlt, rp_errAlt4(i, :))

    grid
    legend("SoT", "RP", "RP (1e-1)");
    xlabel("time [s]");
    ylabel("error [m]");
    switch i
        case 1
            title("X-axis");
        case 2
            title("Y-axis");
        case 3
            title("Z-axis");
    end
end

J4_orerr_fig = figure2('Name', '4th joint orientation errors comparison');
for i = 4:6
    sp_orerr4 = subplot(3, 1, i-3);
    hold on
    
    plot(t_sot, rad2deg(rem(sot_err4(i, :), 2*pi)))
    plot(t_rp, rad2deg(rem(rp_err4(i, :), 2*pi)))
    plot(t_rpAlt, rad2deg(rem(rp_errAlt4(i, :), 2*pi)))

    grid
    legend("SoT", "RP", "RP (1e-1)");
    xlabel("time [s]");
    ylabel("error [deg]");
    switch i
        case 1
            title("X-axis");
        case 2
            title("Y-axis");
        case 3
            title("Z-axis");
    end
end