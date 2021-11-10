clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data
l0 = 0.11; l1 = 0.20; l2 = 0.20; l3 = 0.20; l4 = 0.20; l5 = 0.19; offset = 0.078;

rho = 0.05; % [m]

%% Main trajectory
t_heli = 0:pi/10:3*pi;
x_heli = rho*sin(t_heli);
y_heli = rho*cos(t_heli);
z_heli = 0.6:0.2/30:0.8;

traj_fig = figure2();
grid
hold on
plot3(x_heli, y_heli, z_heli);
for k = 1:size(t_heli,2)
    plot3(x_heli(k), y_heli(k), z_heli(k),'x--','Color','blue')
end
view(50,60)

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

J = @(q)Jac(q,0);
J4 = @(q)Jac(q, 4);

%% Tasks
% Desired EE motion: helicoidal

despos = @(t) [rho*sin(t), rho*cos(t), t, 0, 0, 0];
destwist = @(t) diff(despos(t));

% Desired Joint 4 motion: stable in place
% despos4 = @(t) [0; 0; 0; 0; 0; 0];
% destwist4 = @(t) diff(despos4);

% Oee(t) = [rho*cos(t); rho*sin(t)+l4; l0+l1+l2+l3];%t in terza comp
% OeeHat(t) = AxisToSkew(Oee(t));
% Moee(t) = [[eye(3),     OeeHat(t)];
%            [zeros(3,3), eye(3)   ]];

% Oo4 = [0; 0; l0+l1+l2+l3];
% Oo4Hat = AxisToSkew(Oo4);
% Mo4 = [[eye(3),     Oo4Hat];
%        [zeros(3,3), eye(3)]];

% In global coordinates
% twistEE_O(t) = Moee(t)*twistEE(t);
% twist4_O = Mo4*twist4;

% Complete partial jacobian
J_4 = [J4(q), zeros(size(J4(q),1), size(J(q),2)-size(J4(q),2))];

Rsdes = @(q) RotZ(pi/2);

% Tw(q) = [Twphi1,Twphi2,Twphi3]

% Just use geometric J but compute orientation error with something else 
% (quaternions..) --> look on Sciavicco
TA = [[eye(3),      zeros(3,3)  ];
      [zeros(3,3),  zeros(3,3)]];%Twphi       ]];

Janal = @(q) (TA*J(q));

% jacobs = cat(3,J(q),J_4);

qSotNorm = sot(J(q), ...
               destwist(t)',..., twist4_O],...
               0);

% qSotNorm = sot(J(q), ...
% qSotNorm = sot(Janal(q), ...
%                [twistEE_O(t)],..., twist4_O],...
%                0);

%% Plot solution
f1 = figure2();

for j=[1:size(qSotNorm,1)]
    subplot(size(qSotNorm,1), 1, j)
    plot(squeeze(qSotNorm(j,end,:)))
end

%% Robot design
% TODO: build upon denavit
kuka = rigidBodyTree('Dataformat', 'column');

dhparams = [0  pi/2   0         qSotNorm(1,end,1); 
            0  -pi/2  0         qSotNorm(2,end,1);
            0  -pi/2  l2+l3     qSotNorm(3,end,1);
            0  pi/2   0         qSotNorm(4,end,1);
            0  pi/2   l4+l5     qSotNorm(5,end,1);
            0  -pi/2  0         qSotNorm(6,end,1);
            0  0      offset    qSotNorm(7,end,1)];

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

% f = figure2('Name', 'Using solution');
% % TODO: maybe in a while with a constant check to a flag to repeat?
% % or else, continuously in a thread?
% for i=1:size(qSotNorm,3)
%     gr = show(kuka,qSotNorm(:,end,i));
%     view(-50, 60);
%     drawnow
%     pause(0.1)
% end

%% Integration
ris = cumtrapz(0.01, qSotNorm, 3);

fris = figure2('Name', 'with integration');
for j=1:size(ris,3)
    show(kuka, ris(:,end,3))
    view(-50, 60);
    drawnow
    pause(0.1)
end

hold on
for t=[0:0.1:2*pi]
    plot3(x_heli, y_heli, z_heli);
end

%% Trajectory
% fig_heli = figure2();
% view(-50, 60);
% hold on
% grid
% show(kuka, homeConfiguration(kuka))

fig_traj = figure('Name', 'traj attempt with integration');
for curr_time_step = [1:63]
    timestep = curr_time_step * 0.1;
    xx = double(subs(destwist(t),t,timestep));
    diff_traj(:,curr_time_step) = xx(1:3)';
end

% integration
maybe_traj = cumtrapz(0.1, diff_traj,2);
plot3(maybe_traj(:,1), maybe_traj(:,2), maybe_traj(:,3))