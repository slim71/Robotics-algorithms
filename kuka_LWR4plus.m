%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data
l0 = 11; l1 = 20; l2 = 20; l3 = 20; l4 = 20; l5 = 19; offset = 7.8;

rho = 10; % [cm]

%% Denavit-Hartemberg and Jacobian
%                a   alpha   d       theta    joint_type
denavit = @(q)([[0,  pi/2,   0,      q(1),       "R"] 
                [0,  -pi/2,  0,      q(2),       "R"] 
                [0,  -pi/2,  l2+l3,  q(3),       "R"] 
                [0,  pi/2,   0,      q(4),       "R"] 
                [0,  pi/2,   l4+l5,  q(5),       "R"] 
                [0,  -pi/2,  0,      q(6),       "R"] 
                [0,  0,      offset, q(7),       "R"]]);

Jac = @(q, ind)(JacobFromDH(denavit(q), ind));

J = @(q)Jac(q,0); % J4 = @(q)Jac(q, 4);

Tos = HomX(0, [0,0,l0+l1]);
ForKine = @(q) Tos*DHFKine(denavit(q));

%% Tasks
t_task = 0:pi/10:3*pi;

% Desired EE motion
despose = @(t) [20 - 10*sqrt(3)*sin(t), 20*cos(t), 65 + 10*sin(t),0.5*t,0,sqrt(3)/2*t];

% Make actual derivative of function, avoiding 'diff' problems
destwist = matlabFunction(diff(despose(t)));

traj_fig = figure2('Name', 'Task trajectory');
grid
hold on
first_traj = [];
for k = t_task
    first_traj = [first_traj; despose(k)];
end
plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-')
view(45,45)

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

for k = 1:7
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
curr_pos = @(q) indexAt(ForKine(q),1:3,4);
des_pos = @(t) indexAt(despose(t), 1:3);
pos_error = @(t,q) (des_pos(t) - curr_pos(q)');

curr_or = @(q) MatToQuat(indexAt(ForKine(q),1:3,1:3));
des_or = @(t) EulToQuat(indexAt(despose(t),4:6),'ZYX',false);
or_error = @(t, q) indexAt(QuatProd(des_or(t),QuatInv(curr_or(q))),2:4);

tot_errors = @(t,q) [pos_error(t,q), or_error(t,q)];

%% IK parameters
qinit = [0.5; 0.2; 0.2; pi/3; 0; pi/3; -0.8];

t0 = 0;
tf = 10;

gains_cwe = 0.1 * [1 1 1 1 1 1];

%% Clik kinematics

[t_clik, y_clik] = ode45(@(t,y) clik(t,y,J,destwist), [t0 tf], qinit);
q_clik = y_clik';

ode_fig = figure2('Name', 'Clik with ode45');
for i=[1:size(q_clik, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_clik(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
clik_fig = figure2('Name', 'Clik solution');
for j=[1:size(q_clik,1)]
    subplot(size(q_clik,1), 1, j)
    plot(q_clik(j,:))
end

%% Clik errors
pe_fig = figure2('Name', 'Clik position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_clik,1)]
    pe = pos_error(t_clik(i), q_clik(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_clik,1)]   
    subplot(3, 1, 1);
    plot(t_clik, p1);
    subplot(3, 1, 2);
    plot(t_clik, p2);
    subplot(3, 1, 3);
    plot(t_clik, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'Clik orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_clik,1)]
    oe = or_error(t_clik(i), q_clik(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_clik,1)]   
    subplot(3, 1, 1);
    plot(t_clik, o1);
    subplot(3, 1, 2);
    plot(t_clik, o2);
    subplot(3, 1, 3);
    plot(t_clik, o3);
end

%% Clik with errors
[t_cwe, y_cwe] = ode45(@(t,y) clik_with_error(t,y,J,destwist,tot_errors, gains_cwe),...
                       [t0 tf], qinit);
q_cwe = y_cwe';

ode_fig2 = figure2('Name', 'Cwe with ode45');
for i=[1:size(q_cwe, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_cwe(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
cwe_fig = figure2('Name', 'CWE solution');
for j=[1:size(q_cwe,1)]
    subplot(size(q_cwe,1), 1, j)
    plot(q_cwe(j,:))
end

%% Clik_with_errors errors
pe_fig = figure2('Name', 'CWE position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_cwe,1)]
    pe = pos_error(t_cwe(i), q_cwe(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_cwe,1)]   
    subplot(3, 1, 1);
    plot(t_cwe, p1);
    subplot(3, 1, 2);
    plot(t_cwe, p2);
    subplot(3, 1, 3);
    plot(t_cwe, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'Orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_cwe,1)]
    oe = or_error(t_cwe(i), q_cwe(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_cwe,1)]   
    subplot(3, 1, 1);
    plot(t_cwe, o1);
    subplot(3, 1, 2);
    plot(t_cwe, o2);
    subplot(3, 1, 3);
    plot(t_cwe, o3);
end