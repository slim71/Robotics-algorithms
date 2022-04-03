%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data
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
FK4 = @(q) DHFKine(denavit(q),4);

%% Tasks
t_task = 0:pi/10:3*pi;

% Desired EE motion
despose = @(t) [0.1 * sin(t), 0.2 * cos(t), 0.65, t, t, t]; % + 0.1 * sin(t)
% Make actual derivative of function, avoiding 'diff' problems
destwist = matlabFunction(diff(despose(t)));

% as function of t only for consistency
despose4 = @(t) [0 0 0 0 0 0];
destwist4 = @(t) [0 0 0 0 0 0];

traj_fig = figure2('Name', 'Task trajectory');
title("Trajectory with desone");
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
dhparams = double(denavit(zeros(7,1)));
dhparams = dhparams(:,1:end-1);

bodyNames = {'b1','b2','b3','b4','b5', 'b6','b7'};
parentNames = {'base','b1','b2','b3','b4', 'b5','b6'};
jointNames = {'j1','j2','j3','j4','j5', 'j6','j7'};
jointTypes = {'revolute','revolute','revolute','revolute','revolute', 'revolute','revolute'};
% bodyNames = {'lil','b1','b2','b3','b4','b5', 'b6','b7'};
% parentNames = {'base','lil','b1','b2','b3','b4', 'b5','b6'};
% jointNames = {'j_extra','j1','j2','j3','j4','j5', 'j6','j7'};
% jointTypes = {'fixed', 'revolute','revolute','revolute','revolute','revolute', 'revolute','revolute'};

for k = 1:size(bodyNames,2)
    % Create a rigidBody object with a unique name
    kukaBodies(k) = rigidBody(bodyNames{k});
    % Create a rigidBodyJoint object and give it a unique name
    kukaBodies(k).Joint = rigidBodyJoint(jointNames{k}, jointTypes{k});
    % Use setFixedTransform to specify the body-to-body transformation using DH parameters
%     if k==1
%         setFixedTransform(kukaBodies(k).Joint,Tos);
%         setFixedTransform(kukaBodies(k).Joint,eye(4));
%     else
        setFixedTransform(kukaBodies(k).Joint, dhparams(k,:), 'dh');
%     end
%     if k < 7
%         addVisual(kukaBodies(k),"Mesh", "C:\Program Files\MATLAB\R2022a\toolbox\robotics\robotmanip\exampleRobots\iiwa_description\meshes\iiwa14\visual\link_"+k+".stl")
%     end
    % Attach the body joint to the robot
    addBody(kuka, kukaBodies(k), parentNames{k});
end

showdetails(kuka)

%% Needed functions
curr_pos = @(q) indexAt(ForKine(q),1:3,4);
des_pos = @(t) indexAt(despose(t), 1:3);
pose_error = @(t,q) (des_pos(t) - curr_pos(q)');

curr_or = @(q) MatToQuat(indexAt(ForKine(q),1:3,1:3));
des_or = @(t) EulToQuat(indexAt(despose(t),4:6),'ZYX',false);
or_error = @(t, q) indexAt(QuatProd(des_or(t),QuatInv(curr_or(q))),2:4);

tot_errors = @(t,q) [pose_error(t,q), or_error(t,q)];

%% IK parameters
qinit = [-pi/2, pi/4, 0, pi/3,0,0,0];
% qinit = [-pi/3, pi/4, pi/4, 2*pi/3, -pi/6, 3*pi/4, pi/3];

t0 = 0;
tf = 7;

%% SOT 1 task
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
for k=[1:size(q1_sot, 2)]
    fk = ForKine(q1_sot(:, k));
    positions = fk(1:3, 4);
    plot3(positions(1), positions(2), positions(3), 'o');
end

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

%% SoT
[t_sot, y_sot] = ode15s(@(t,y) sot(t, y, {J, J4}, {task1 destwist4}, 10^(-1)),...
                       [t0 tf], qinit, options);

q_sot = y_sot';

sot_fig = figure2('Name', 'SoT with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
for i=[1:size(q_sot, 2)]
    plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot(:, i));
    view(45, 45)
    drawnow
    pause(0.1)
    hold off
end

% sot_fig2 = figure2('Name', 'SoT with ode15s');
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

%% Reverse Priority
%{
[t_rp, y_rp] = ode45(@(t,y) rp(t,y,{J, J4},{task1 destwist4}, 10^(-1)),...
                       [t0 tf], qinit);

q_rp = y_rp';

rp_fig = figure2('Name', 'RP with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
plot3(first_traj(:, 1), first_traj(:, 2), first_traj(:, 3), '-', 'Color', 'black');
hold on
show(kuka,qinit');
grid
for k=[1:size(q_sotone,2)]
    fk = ForKine(q_sotone(:,k));
    posi = fk(1:3,4);
    plot3(posi(1), posi(2),posi(3),'o');
end

rp_fig2 = figure2('Name', 'RP with ode15s');
title('$$ \lambda = 10^{-1} $$', 'interpreter', 'latex');
for i=[1:size(q_rp, 2)]
    plot3(first_traj(:, 1), first_traj(:, 2),first_traj(:, 3), '-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp(:, i));
    view(45, 45)
    drawnow
    pause(0.1)
    hold off
end

%}