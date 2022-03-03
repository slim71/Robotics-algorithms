%% Preparations
clc; clear all; close all;

syms t
syms q [7 1]
% TODO: assumptions?

%% Numerical Data
l0 = 11; l1 = 20; l2 = 20; l3 = 20; l4 = 20; l5 = 19; offset = 7.8;

% rho = 10; % [cm]

%% Denavit-Hartemberg and Jacobian
%                a   alpha   d       theta    joint_type
denavit = @(q)([[0,  pi/2,   0,      q(1),       "R"] 
                [0,  -pi/2,  0,      q(2),       "R"] 
                [0,  -pi/2,  l2+l3,  q(3),       "R"] 
                [0,  pi/2,   0,      q(4),       "R"] 
                [0,  pi/2,   l4+l5,  q(5),       "R"] 
                [0,  -pi/2,  0,      q(6),       "R"] 
                [0,  0,      offset, q(7),       "R"]]);

J = @(q) DHJacob0(denavit(q), 0);
J4 = @(q) DHJacob0(denavit(q), 4);

Tos = HomX(0, [0,0,l0+l1]);

ForKine = @(q) Tos*DHFKine(denavit(q));
FK4 = @(q) Tos*DHFKine(denavit(q),4);

%% Tasks
t_task = 0:pi/10:3*pi;

% Desired EE motion
despose = @(t) [20 - 10*sqrt(3)*sin(t), 20*cos(t), 65 + 10*sin(t),0.5*t,0,sqrt(3)/2*t];
% Make actual derivative of function, avoiding 'diff' problems
destwist = matlabFunction(diff(despose(t)));

% as function of t only for consistency
despose4 = @(t) [0 0 0 0 0 0];
destwist4 = @(t) [0 0 0 0 0 0];

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

bodyNames = {'lil','b1','b2','b3','b4','b5', 'b6','b7'};
parentNames = {'base','lil','b1','b2','b3','b4', 'b5','b6'};
jointNames = {'j_extra','j1','j2','j3','j4','j5', 'j6','j7'};
jointTypes = {'fixed', 'revolute','revolute','revolute','revolute','revolute', 'revolute','revolute'};

for k = 1:size(bodyNames,2)
    % Create a rigidBody object with a unique name
    kukaBodies(k) = rigidBody(bodyNames{k});
    % Create a rigidBodyJoint object and give it a unique name
    kukaBodies(k).Joint = rigidBodyJoint(jointNames{k}, jointTypes{k});
    % Use setFixedTransform to specify the body-to-body transformation using DH parameters
    if k==1
        setFixedTransform(kukaBodies(k).Joint,Tos);
    else
        setFixedTransform(kukaBodies(k).Joint, dhparams(k-1,:), 'dh');
    end
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
% qinit = [0.5; 0.2; 0.2; pi/3; 0; pi/3; -0.8];
% qinit = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
% {-1.18276, 0.899385, 0.789808, 1.97911, -0.32927, 2.44162, -1.48671}
qinit = [-pi/3, pi/4, pi/4, 2*pi/3, -pi/6, 3*pi/4, pi/3];

t0 = 0;
tf = 7;

%% Clik kinematics

[t_clik, y_clik] = ode15s(@(t,y) clik(t,y,J,destwist), [t0 tf], qinit);
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
gains_cwe = 50 * [1 1 1 1 1 1];

[t_cwe, y_cwe] = ode15s(@(t,y) clik_with_error(t,y,J,destwist,tot_errors, gains_cwe),...
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

or_fig = figure2('Name', 'CWE orientation error');
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
%% Clik with secondary objective
gains_cso = 50 * [1 1 1 1 1 1];

k0 = eye(7);
mins = [-pi, -pi, -pi, -pi, -pi, -pi, -pi];
maxs = [pi, pi, pi, pi, pi, pi, pi];

H = @(q) JointCenterRange(q, mins, maxs);

[t_cso, y_cso] = ode45(@(t,y) clik_secobj(t,y,J,destwist,tot_errors, gains_cso, H),...
                       [t0 tf], qinit);
q_cso = y_cso';

ode_fig2 = figure2('Name', 'Cso with ode45');
for i=[1:size(q_cso, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_cso(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
cwe_fig = figure2('Name', 'CSO solution');
for j=[1:size(q_cso,1)]
    subplot(size(q_cso,1), 1, j)
    plot(q_cso(j,:))
end

%% Clik_secobj errors
pe_fig = figure2('Name', 'CSO position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_cso,1)]
    pe = pos_error(t_cso(i), q_cso(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_cso,1)]   
    subplot(3, 1, 1);
    plot(t_cso, p1);
    subplot(3, 1, 2);
    plot(t_cso, p2);
    subplot(3, 1, 3);
    plot(t_cso, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'CSO Orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_cso,1)]
    oe = or_error(t_cso(i), q_cso(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_cso,1)]   
    subplot(3, 1, 1);
    plot(t_cso, o1);
    subplot(3, 1, 2);
    plot(t_cso, o2);
    subplot(3, 1, 3);
    plot(t_cso, o3);
end

%% Clik with q halfway
gains_half = 50 * [1 1 1 1 1 1];
mins = [-pi, -pi, -pi, -pi, -pi, -pi, -pi];
maxs = [pi, pi, pi, pi, pi, pi, pi];

[t_half, y_half] = ode45(@(t,y) clik_secobj(t,y,J,destwist,tot_errors, gains_half, mins, maxs),...
                       [t0 tf], qinit);
q_half = y_half';

ode_half = figure2('Name', 'Halfq with ode');
for i=[1:size(q_half, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_half(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
half_fig = figure2('Name', 'CSO solution');
for j=[1:size(q_half,1)]
    subplot(size(q_half,1), 1, j)
    plot(q_half(j,:))
end

%% Clik_halfq errors
pe_fig = figure2('Name', 'CSO position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_half,1)]
    pe = pos_error(t_half(i), q_half(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_half,1)]   
    subplot(3, 1, 1);
    plot(t_half, p1);
    subplot(3, 1, 2);
    plot(t_half, p2);
    subplot(3, 1, 3);
    plot(t_half, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'Halfq Orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_half,1)]
    oe = or_error(t_half(i), q_half(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_half,1)]   
    subplot(3, 1, 1);
    plot(t_half, o1);
    subplot(3, 1, 2);
    plot(t_half, o2);
    subplot(3, 1, 3);
    plot(t_half, o3);
end

%% Manipolability index TODO later
% manIndex = @(q) sqrt(det(J(q) * J(q)'));
% 
% H2 = matlabFunction(diff(manIndex(q)));dX

%% SoT
% options = odeset('RelTol', 1e-2, 'AbsTol', 1e-3, 'InitialStep', 1e-2);

[t_sot6, y_sot6] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-6)),...
                       [t0 tf], qinit);
[t_sot5, y_sot5] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-5)),...
                       [t0 tf], qinit);
[t_sot4, y_sot4] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-4)),...
                       [t0 tf], qinit);
[t_sot3, y_sot3] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-3)),...
                       [t0 tf], qinit);
[t_sot2, y_sot2] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-2)),...
                       [t0 tf], qinit);
[t_sot1, y_sot1] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(-1)),...
                       [t0 tf], qinit);
[t_sot0, y_sot0] = ode45(@(t,y) sot(t,y,{J, J4},{destwist destwist4}, 10^(0)),...
                       [t0 tf], qinit);

q_sot6 = y_sot6';

sot_fig6 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-6} $$','interpreter','latex');
for i=[1:size(q_sot6, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot6(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot5 = y_sot5';

sot_fig5 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-5} $$','interpreter','latex');
for i=[1:size(q_sot5, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot5(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot4 = y_sot4';

sot_fig4 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-4} $$','interpreter','latex');
for i=[1:size(q_sot4, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot4(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot3 = y_sot3';

sot_fig3 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-3} $$','interpreter','latex');
for i=[1:size(q_sot3, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot3(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot2 = y_sot2';

sot_fig2 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-2} $$','interpreter','latex');
for i=[1:size(q_sot2, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot2(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot1 = y_sot1';

sot_fig1 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{-1} $$','interpreter','latex');
for i=[1:size(q_sot1, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot1(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_sot0 = y_sot0';

sot_fig0 = figure2('Name', 'SoT with ode45');
title('$$ \lambda = 10^{0} $$','interpreter','latex');
for i=[1:size(q_sot0, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_sot0(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

%% Reverse Priority
[t_rp, y_rp] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-6)),...
                       [t0 tf], qinit);
[t_rp2, y_rp2] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-5)),...
                       [t0 tf], qinit);
[t_rp3, y_rp3] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-4)),...
                       [t0 tf], qinit);
[t_rp4, y_rp4] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-3)),...
                       [t0 tf], qinit);
[t_rp5, y_rp5] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-2)),...
                       [t0 tf], qinit);
[t_rp6, y_rp6] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(-1)),...
                       [t0 tf], qinit);
[t_rp7, y_rp7] = ode45(@(t,y) rp(t,y,{J, J4},{destwist destwist4}, 10^(0)),...
                       [t0 tf], qinit);

q_rp = y_rp';

rp_fig = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-6} $$','interpreter','latex');
for i=[1:size(q_rp, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp2 = y_rp2';

rp_fig2 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-5} $$','interpreter','latex');
for i=[1:size(q_rp2, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp2(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp3 = y_rp3';

rp_fig3 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-4} $$','interpreter','latex');
for i=[1:size(q_rp3, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp3(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp4 = y_rp4';

rp_fig4 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-3} $$','interpreter','latex');
for i=[1:size(q_rp4, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp4(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp5 = y_rp5';

rp_fig5 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-2} $$','interpreter','latex');
for i=[1:size(q_rp5, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp5(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp6 = y_rp6';

rp_fig6 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{-1} $$','interpreter','latex');
for i=[1:size(q_rp6, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp6(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end

q_rp7 = y_rp7';

rp_fig7 = figure2('Name', 'RP with ode45');
title('$$ \lambda = 10^{0} $$','interpreter','latex');
for i=[1:size(q_rp7, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_rp7(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end