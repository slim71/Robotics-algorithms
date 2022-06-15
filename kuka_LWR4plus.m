%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]

error_titles = [
    "X-axis position";
    "Y-axis position";
    "Z-axis position";
    "X-axis orientation";
    "Y-axis orientation";
    "Z-axis orientation"
    ];
%% Numerical Data [cm]
% Robot measurements
l0 = 0.11; l1 = 0.20; l2 = 0.20; l3 = 0.20; l4 = 0.20; l5 = 0.19;

%Circle radius
r = 0.2;

% Task height from ground
h = 0.65;

% Dampening factor
lambda = 1e-3;

%% IK parameters
qinit = [-pi/2, pi/4, 0, pi/3,0,0,0];

t0 = 0;
tf = 7;
t_circle = 0:pi/10:3*pi;

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

%% End-effector functions
% Using matlabFunction(diff()) to make actual derivative of function, 
% avoiding 'diff' classic problems.

% Desired EE position
EE_des_pos = @(t) [r * sin(t), r * cos(t), h];

% Tanget and normal versor to the EE task
% (and subsequent perpendicular axis)
tan_vect = matlabFunction(diff(EE_des_pos(t)));
vect_module = @(t) sqrt(indexAt(tan_vect(t), 1)^2 + ...
                       indexAt(tan_vect(t),2)^2 + ...
                       indexAt(tan_vect(t),3)^2);
tan_vers = @(t) tan_vect(t)/vect_module(t);
norm_vect = matlabFunction(diff(tan_vers(t)));
norm_module = @(t) sqrt(indexAt(norm_vect(t), 1)^2 + ...
                       indexAt(norm_vect(t),2)^2 + ...
                       indexAt(norm_vect(t),3)^2);
norm_vers = @(t) norm_vect(t)/norm_module(t);
z_vers = @(t) cross(tan_vers(t), norm_vers(t));

% Desired EE orientation
EE_Rdes = @(t) [tan_vers(t)', norm_vers(t)', z_vers(t)'];
EE_des_or = @(t) MatToRPY(EE_Rdes(t));

% Desired EE pose
EE_despose = @(t) [EE_des_pos(t), EE_des_or(t)];

% Desired EE twist
EE_v_des = matlabFunction(diff(EE_des_pos(t))); % linear velocity
EE_Rdesdot = matlabFunction(diff(EE_Rdes(t)));
EE_w_des = @(t) HatInv(EE_Rdesdot(t)*EE_Rdes(t)'); % angular velocity
EE_destwist = @(t) [EE_v_des(t), EE_w_des(t)];

% EE position error
EE_curr_pos = @(q) indexAt(ForKine(q), 1:3, 4)';
EE_pos_error = @(t, q) (EE_des_pos(t) - EE_curr_pos(q));

% EE orientation error
EE_curr_or_quat = @(q) MatToQuat(indexAt(ForKine(q), 1:3, 1:3));
EE_des_or_quat = @(t) EulToQuat(indexAt(EE_despose(t), 4:6), 'ZYX', false);
EE_or_error = @(t, q) indexAt(QuatProd(EE_des_or_quat(t), QuatInv(EE_curr_or_quat(q))), 2:4);

% EE pose error
EE_pose_errors = @(t, q) [EE_pos_error(t, q), EE_or_error(t, q)];

%% 4th joint functions
% Some functions have t and q as arguments only for consistency

% 4th joint motion
J4_des_pos = @(t) indexAt(FK4(qinit), 1:3, 4)';
J4_des_or = @(t) MatToEulZYX(indexAt(FK4(qinit), 1:3, 1:3));
J4_despose = @(t) [J4_des_pos(t), J4_des_or(t)];

% J4_v_des = @(t) [0, 0, 0];
% J4_w_des = @(t) [0, 0, 0]; % even using Rdesdot, it would be 0s
% J4_destwist4 = @(t) [J4_v_des, J4_w_des];
J4_destwist4 = @(t, q) [0 0 0 0 0 0];

J4_curr_pos = @(q) indexAt(FK4(q), 1:3, 4)';
% J4_des_pos = @(t) indexAt(J4_despose(t), 1:3);
J4_pos_error = @(t, q) (J4_des_pos(t) - J4_curr_pos(q));

J4_curr_or_quat = @(q) MatToQuat(indexAt(FK4(q), 1:3, 1:3));
J4_des_or_quat = @(t) EulToQuat(indexAt(J4_despose(t), 4:6), 'ZYX', false);
J4_or_error = @(t, q) indexAt(QuatProd(J4_des_or_quat(t), QuatInv(J4_curr_or_quat(q))), 2:4);

J4_pose_errors = @(t, q) [J4_pos_error(t, q), J4_or_error(t, q)];

%% Tasks

[des_circle_x, des_circle_y, des_circle_z] = circleFromFun(EE_des_pos, t_circle);

% Graphical representation of the circle describing the EE main task
traj_fig = figure2('Name', 'Task trajectory');
plot3(des_circle_x, des_circle_y, des_circle_z, 'k-');
title("Main task trajectory");
grid on
hold on
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-2*r, 2*r]); ylim([-2*r, 2*r]); zlim([0, 2*h]);
view(45, 45);

% Draw an example of the desired EE orientation
init_p = indexAt(EE_despose(0),1:3);
start_x = tan_vers(0)' * 0.1;
start_y = norm_vers(0)' * 0.1;
start_z = z_vers(0)' * 0.2;

EEOrientation(init_p, [start_x, start_y, start_z]);

% Mathematical formulas for the two tasks
K = 50 * eye(6);
K2 = 50 * eye(6);
task1 = @(t, q) (K * EE_pose_errors(t, q)')';
task2 = @(t, q) (K2 * J4_pose_errors(t, q)')';

%% Robot design
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

%% SoT
[t_sot, y_sot] = ode15s(@(t,y) sot(t, y, {J, J4}, {task1, task2}, lambda),...
                       [t0 tf], qinit);

q_sot = y_sot';

sot_fig = figure2('Name', 'SoT with ode15s');
plot3(des_circle_x, des_circle_y, des_circle_z, 'k-');
title('SoT, $$ \lambda = 10^{-3} $$', 'interpreter', 'latex');
hold on
show(kuka, qinit');
grid
for k = 1:size(q_sot, 2)
    fk = ForKine(q_sot(:, k));
    sot_positions = fk(1:3, 4);
    plot3(sot_positions(1), sot_positions(2), sot_positions(3), 'o');
end

% Animated plot
sot2_anim = figure2('Name', 'SoT with ode15s');
for i = 1:size(q_sot, 2)
    % Circle
    plot3(des_circle_x, des_circle_y, des_circle_z, 'k-');
    title('SoT, $$ \lambda = 10^{-3} $$', 'interpreter', 'latex');
    hold on
    grid
    % Robot and EE position
    s_sot = show(kuka, q_sot(:, i), 'Visuals', 'off');
    patch_list = findall(s_sot, "Type", "Patch");
    for x = 1:size(patch_list, 1)
        patch_list(x).Visible = 'off';
    end
    % EE orientation
    sot_R = indexAt(ForKine(q_sot(:, i)), 1:3, 1:3);
    EEOrientation(indexAt(ForKine(q_sot(:, i)), 1:3, 4), sot_R*0.025);
    view(45, 45)
    drawnow
    pause(0.1)
    hold off
end

%% SoT errors

sot_err = zeros(6, size(q_sot, 2));
sot_err4 = zeros(6, size(q_sot, 2));
for k = 1:size(q_sot, 2)
    sot_err(:, k) = EE_pose_errors(t_sot(k), q_sot(:, k));
    sot_err4(:, k) = J4_pose_errors(t_sot(k), q_sot(:, k));
end

soterr_fig = figure2('Name', 'SoT errors');
sgtitle("SoT errors");
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
    title(error_titles(i));
end

%% Reverse Priority
% not much different from 10^-3 on
[t_rp, y_rp] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, task2}, 10^(-3), 0),...
                       [t0 tf], qinit);

q_rp = y_rp';

rp_fig = figure2('Name', 'RP with ode15s');
plot3(des_circle_x, des_circle_y, des_circle_z, 'k-');
title('RP, $$ \lambda = 10^{-3} $$', 'interpreter', 'latex');
hold on
s = show(kuka, qinit');
grid
for k = 1:size(q_rp, 2)
    fk = ForKine(q_rp(:, k));
    rp_position = fk(1:3, 4);
    plot3(rp_position(1), rp_position(2), rp_position(3), 'o');
end

% Animated plot
rp_fig2 = figure2('Name', 'RP with ode15s');
for i= 1:size(q_rp, 2)
    plot3(des_circle_x, des_circle_y, des_circle_z, 'k-');
    title('RP, $$ \lambda = 10^{-3} $$', 'interpreter', 'latex');
    hold on
    grid on
    s_rp = show(kuka, q_rp(:, i));
    patch_list = findall(s_rp, "Type", "Patch");
    for x = 1:size(patch_list, 1)
        patch_list(x).Visible = 'off';
    end
    % EE orientation
    rp_R = indexAt(ForKine(q_rp(:, i)), 1:3, 1:3);
    EEOrientation(indexAt(ForKine(q_rp(:, i)), 1:3, 4), rp_R*0.025);
    view(45, 45)
    drawnow
    pause(0.1)
    hold off
end

%% RP errors

rp_err = zeros(6, size(q_rp, 2));
rp_err4 = zeros(6, size(q_rp, 2));
for k = 1:size(q_rp, 2)
    rp_err(:, k) = EE_pose_errors(t_rp(k), q_rp(:, k));
    rp_err4(:, k) = J4_pose_errors(t_rp(k), q_rp(:, k));
end

rperr_fig = figure2('Name', 'RP errors');
sgtitle("RP errors");
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
    title(error_titles(i));
end

%% RP with Rank-one update
[t_rprou, y_rprou] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, task2}, 10^(-3), 1),...
                       [t0 tf], qinit);

q_rprou = y_rprou';

rprou_err = zeros(6, size(q_rprou, 2));
rprou_err4 = zeros(6, size(q_rprou, 2));
for k = 1:size(q_rprou, 2)
    rprou_err(:, k) = EE_pose_errors(t_rprou(k), q_rprou(:, k));
    rprou_err4(:, k) = J4_pose_errors(t_rprou(k), q_rprou(:, k));
end

rprou_poserr_fig = figure2('Name', 'RP (w/ rank-one update) EE errors');
sgtitle("RP (w/ rank-one update) EE errors");
for i = 1:6
    sp_rpo1 = subplot(6, 1, i);
    hold on

    if i <= 3 % position errors
        plot(t_rp, rp_err(i, :))
        plot(t_rprou, rprou_err(i, :))

        ylabel("error [m]");

    else % orientation error
        plot(t_rp, rad2deg(rem(rp_err(i, :), 2*pi)))
        plot(t_rprou, rad2deg(rem(rprou_err(i, :), 2*pi)))

        ylabel("error [deg]");

    end

    grid
    legend("Normal RP", "With rank-one update");
    xlabel("time [s]");
    title(error_titles(i));
end

rprou_orerr_fig = figure2('Name', 'RP (w/ rank-one update) J4 errors');
sgtitle("RP (w/ rank-one update) J4 errors");
for i = 1:6
    sp_rpo1 = subplot(6, 1, i);
    hold on

    if i <= 3 % position errors
        plot(t_rp, rp_err4(i, :))
        plot(t_rprou, rprou_err4(i, :))

        ylabel("error [m]");

    else % orientation error
        plot(t_rp, rad2deg(rem(rp_err4(i, :), 2*pi)))
        plot(t_rprou, rad2deg(rem(rprou_err4(i, :), 2*pi)))

        ylabel("error [deg]");

    end

    grid
    legend("Normal RP", "With rank-one update");
    xlabel("time [s]");
    title(error_titles(i));
end

%% E-E error comparison

poserr_fig = figure2('Name', 'EE position errors comparison');
sgtitle("EE position errors")
for i = 1:3
    subplot(3, 1, i)
    hold on

    plot(t_sot, sot_err(i, :))
    plot(t_rp, rp_err(i, :))
    plot(t_rprou, rprou_err(i, :))

    grid
    legend("SoT", "RP", "RP-ROU");
    xlabel("time [s]");
    ylabel("error [m]");
    title(error_titles(i));
end

% for i=1:3
%     subplot(3,1,i);
%     xlim([2.5, 3.5])
% end

orerr_fig = figure2('Name', 'EE orientation errors comparison');
sgtitle("EE orientation errors")
for i = 4:6
    sp_orerr = subplot(3, 1, i-3);
    hold on

    plot(t_sot, rad2deg(rem(sot_err(i, :), 2*pi)));
    plot(t_rp, rad2deg(rem(rp_err(i, :), 2*pi)));
    plot(t_rprou, rad2deg(rem(rprou_err(i, :), 2*pi)));

    grid
    legend("SoT", "RP", "RP-ROU");
    xlabel("time [s]");
    ylabel("error [deg]");
    title(error_titles(i));
end

% for i=1:3
%     subplot(3,1,i);
%     xlim([1.5, 2.5])
% end

%% 4th joint error comparison

J4_poserr_fig = figure2('Name', '4th joint position errors comparison');
sgtitle("4th joint position errors");
for i = 1:3
    subplot(3, 1, i)
    hold on

    plot(t_sot, sot_err4(i, :))
    plot(t_rp, rp_err4(i, :))
    plot(t_rprou, rprou_err4(i, :))

    grid
    legend("SoT", "RP", "RP-ROU");
    xlabel("time [s]");
    ylabel("error [m]");
    title(error_titles(i));
end

% for i=1:3
%     subplot(3,1,i);
%     xlim([1.5, 2])
% end

J4_orerr_fig = figure2('Name', '4th joint orientation errors comparison');
sgtitle("4th joint orientation errors");
for i = 4:6
    sp_orerr4 = subplot(3, 1, i-3);
    hold on
    
    plot(t_sot, rad2deg(rem(sot_err4(i, :), 2*pi)))
    plot(t_rp, rad2deg(rem(rp_err4(i, :), 2*pi)))
    plot(t_rprou, rad2deg(rem(rprou_err4(i, :), 2*pi)))

    grid
    legend("SoT", "RP", "RP-ROU");
    xlabel("time [s]");
    ylabel("error [deg]");
    title(error_titles(i));
end

% for i=1:3
%     subplot(3,1,i);
%     xlim([1.5, 2.5])
% end
%% Usefule notes
time_loops = 10;
sot_times = zeros(1, time_loops);
rp_times = zeros(1, time_loops);
rprou_times = zeros(1, time_loops);

for count = 1:time_loops
    tic
    [~, ~] = ode15s(@(t,y) sot(t, y, {J, J4}, {task1, task2}, lambda),...
                           [t0 tf], qinit);
    sot_times(count) = toc;
    tic
    [~, ~] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, task2}, 10^(-3), 0),...
                           [t0 tf], qinit);
    rp_times(count) = toc;
    tic
    [~, ~] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, task2}, 10^(-3), 1),...
                           [t0 tf], qinit);
    rprou_times(count) = toc;


end

methods = ["SoT"; "RP"; "RP-ROU"];
comp_time = [mean(sot_times); mean(rp_times); mean(rprou_times)];
mean_sot_error = mean(sot_err, 2);
mean_rp_error = mean(rp_err, 2);
mean_rprou_error = mean(rprou_err, 2);
max_sot_error = max(sot_err, [], 2);
max_rp_error = max(rp_err, [], 2);
max_rprou_error = max(rprou_err, [], 2);

table(comp_time, ...
    [mean_sot_error(1:3)'; mean_rp_error(1:3)'; mean_rprou_error(1:3)'], ...
    [rad2deg(rem(mean_sot_error(4:6)', 2*pi)); 
        rad2deg(rem(mean_rp_error(4:6)', 2*pi)); 
        rad2deg(rem(mean_rprou_error(4:6)', 2*pi))], ...
    [max_sot_error(1:3)'; max_rp_error(1:3)'; max_rprou_error(1:3)'], ...
    [rad2deg(rem(max_sot_error(4:6)', 2*pi));
        rad2deg(rem(max_rp_error(4:6)', 2*pi)); 
        rad2deg(rem(max_rprou_error(4:6)', 2*pi))], ...
     'VariableNames', ...
     {'Computational time', 'Mean position errors (x y z)', ...
     'Mean orientation errors (x y z)', 'Max position errors (x y z)', ...
     'Max orientation errors (x y z)'}, ...
     'RowNames', methods)
   