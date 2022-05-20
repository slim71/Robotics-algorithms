%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data [cm]
l0 = 0.11; l1 = 0.20; l2 = 0.20; l3 = 0.20; l4 = 0.20; l5 = 0.19;

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
EE_des_pos = @(t) [0.2 * sin(t), 0.2 * cos(t), 0.65];

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

% Graphical representation of the circle describing the EE main task
traj_fig = figure2('Name', 'Task trajectory');
title("Main task trajectory");
grid
hold on
first_traj = [];
for k = t_circle
    first_traj = [first_traj; EE_despose(k)];
end
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-')
xlabel('x');ylabel('y');zlabel('z');
view(45, 45)

% Draw an example of the desired EE orientation
init_p = indexAt(EE_despose(0),1:3);
x_dir = tan_vers(k)';
y_dir = norm_vers(k)';
z_dir = z_vers(k)';

quiv_x = quiver3(init_p(1), init_p(2), init_p(3), x_dir(1), x_dir(2), x_dir(3));
quiv_x.Color = 'red';
quiv_y = quiver3(init_p(1), init_p(2), init_p(3), y_dir(1), y_dir(2), y_dir(3));
quiv_y.Color = 'green';
quiv_z = quiver3(init_p(1), init_p(2), init_p(3), z_dir(1), z_dir(2), z_dir(3));
quiv_z.Color = 'blue';

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
[t_sot, y_sot] = ode15s(@(t,y) sot(t, y, {J, J4}, {task1, task2}, 10^(-3)),...
                       [t0 tf], qinit);

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
    sot_err(:, k) = EE_pose_errors(t_sot(k), q_sot(:, k));
    sot_err4(:, k) = J4_pose_errors(t_sot(k), q_sot(:, k));
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

%% Reverse Priority

% not much different from 10^-3 on
[t_rp, y_rp] = ode45(@(t,y) rp(t, y, {J, J4}, {task1, task2}, 10^(-3)),...
                       [t0 tf], qinit);

q_rp = y_rp';

rp_fig = figure2('Name', 'RP with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka, qinit');
grid
for k = 1:size(q_rp, 2)
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
    rp_err(:, k) = EE_pose_errors(t_rp(k), q_rp(:, k));
    rp_err4(:, k) = J4_pose_errors(t_rp(k), q_rp(:, k));
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