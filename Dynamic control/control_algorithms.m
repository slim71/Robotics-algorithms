%% Preparations
clc; close all; clear all

qlabels = ["$q_{1} [rad]$"; "$q_{2} [rad]$"; "$q_{3} [rad]$"; "$q_{4} [rad]$"; "$q_{5} [rad]$"; "$q_{6} [rad]$"];
qdlabels = ["$\dot{q_{1}} [rad/s]$"; "$\dot{q_{2}} [rad/s]$"; "$\dot{q_{3}} [rad/s]$"; 
            "$\dot{q_{4}} [rad/s]$"; "$\dot{q_{5}} [rad/s]$"; "$\dot{q_{6}} [rad/s]$"];
qddlabels = ["$\ddot{q_{1}} [rad/s^2]$"; "$\ddot{q_{2}} [rad/s^2]$"; "$\ddot{q_{3}} [rad/s^2]$"; 
             "$\ddot{q_{4}} [rad/s^2]$"; "$\ddot{q_{5}} [rad/s^2]$"; "$\ddot{q_{6}} [rad/s^2]$"];
eeaxislabels = ["x [mm]"; "y [mm]"; "z [mm]"; "x [deg]"; "y [deg]"; "z [deg]"];

%% Importing needed robot model

abbirb_urdf = importrobot("abbIrb1600.urdf");

w_mesh = figure2('Name', 'IRB1660 model with mesh');
show(abbirb_urdf, 'visuals', 'on', 'collision', 'off');
xlim([-0.5, 1]);
ylim([-0.5, 0.5]);
zlim([0, 1.1]);
ylabel("Y [m]");
xlabel("X [m]");
zlabel("Z [m]");

wo_mesh = figure2('Name', 'IRB1660 model without mesh');
show(abbirb_urdf, 'visuals', 'off', 'collision', 'off');
xlim([-0.5, 1]);
ylim([-0.5, 0.5]);
zlim([0, 1.1]);

%% Robot dimensions and DH parameters

% Masses [g]
% Links 3 through six are treated as a single body
m1 = 104393.43;
m2 = 20452.12;
m3 = 50588.7;
% Since the total mass of the robot is 250kg from the manufacturer's
% datasheets, we distribute the rest of the mass on the last links
m456 = (250e3 - (m1 + m2 + m3)) / 3;

% Center of masses [mm]
x_cm1 = 51.98;
x_cm2 = 201.93;
x_cm3 = -11.12;
y_cm1 = -343.26;
y_cm2 = -0.15;
y_cm3 = 15.28;
z_cm1 = -11.92;
z_cm2 = -182.8;
z_cm3 = 96.81;

% Lenghts [mm]
d1 = 486.5;
d4 = 600;
d6 = 65;
a1 = 150;
a2 = 475;

% Angles [rad]
alpha1 = -pi/2;
alpha3 = -pi/2;
alpha4 = pi/2;
alpha5 = pi/2;
theta2_0 = -pi/2;
theta5_0 = pi;

%% Denavit-Hartenberg
% 6 joints, 5 links + 1 E-E

% Links
L1 = Link('d', d1, 'a', a1, 'alpha', alpha1, 'offset', 0,        'revolute');
L2 = Link('d', 0,  'a', a2, 'alpha', 0,      'offset', theta2_0, 'revolute');
L3 = Link('d', 0,  'a', 0,  'alpha', alpha3, 'offset', 0,        'revolute');
L4 = Link('d', d4, 'a', 0,  'alpha', alpha4, 'offset', 0,        'revolute');
L5 = Link('d', 0,  'a', 0,  'alpha', alpha5, 'offset', theta5_0, 'revolute');
L6 = Link('d', d6, 'a', 0,  'alpha', 0,      'offset', 0,        'revolute');

L1_unc = Link('d', d1, 'a', a1, 'alpha', alpha1, 'offset', 0,        'revolute');
L2_unc = Link('d', 0,  'a', a2, 'alpha', 0,      'offset', theta2_0, 'revolute');
L3_unc = Link('d', 0,  'a', 0,  'alpha', alpha3, 'offset', 0,        'revolute');
L4_unc = Link('d', d4, 'a', 0,  'alpha', alpha4, 'offset', 0,        'revolute');
L5_unc = Link('d', 0,  'a', 0,  'alpha', alpha5, 'offset', theta5_0, 'revolute');
L6_unc = Link('d', d6, 'a', 0,  'alpha', 0,      'offset', 0,        'revolute');

%% Robot models

% Real robot
abbirb = SerialLink([L1, L2, L3, L4, L5, L6], 'name', 'ABBIRB1600', 'plotopt', {'notiles'});

% Uncertan model of the robot (used in adaptive controls) 
abbirb_unc = SerialLink([L1, L2, L3, L4, L5, L6], 'name', 'ABBIRB1600unc', 'plotopt', {'notiles'});

%% Set real link masses, lengths and inertia
% ref. system: z up, x forward

% Masses
abbirb.links(1).m = m1;
abbirb.links(2).m = m2;
abbirb.links(3).m = m3;
abbirb.links(4).m = m456;
abbirb.links(5).m = m456;
abbirb.links(6).m = m456;

% Centers of gravity; links 3 through 6 are considered as together
abbirb.links(1).r = [x_cm1, y_cm1, z_cm1];
abbirb.links(2).r = [x_cm2, y_cm2, z_cm2];
abbirb.links(3).r = [x_cm3, y_cm3, z_cm3];
% We suppose the last links are all "together"
abbirb.links(4).r = abbirb.links(3).r;
abbirb.links(5).r = abbirb.links(3).r;
abbirb.links(6).r = abbirb.links(3).r;

% Inertia matrices [g * mm^2]
abbirb.links(1).I = [[1077815.9,     1939.4,         24207.0];
                     [1939.4,        1010823.7,     -13131.7];
                     [24207.0,      -13131.7,        225407.6]];
abbirb.links(2).I = [[688775026.5,  -2153417.8,     -83961.2];
                     [-2153417.8,    52950060.9,    -48598551.6];
                     [-83961.2,     -48598551.6,     703024081.4]];
abbrib.links(3).I = [[249444.6,      0.0,            0.0];
                     [0.0,           305986.7,      -0.1];
                     [0.0,          -0.1,           388391.5]];
abbirb.links(4).I = abbrib.links(3).I;
abbirb.links(5).I = abbrib.links(3).I;
abbirb.links(6).I = abbrib.links(3).I;

%% Uncertain model parameters

% Percentage of uncertainty on parameters value
unc = 10; % 10%

% Masses
abbirb_unc.links(1).m = m1 * (1 + unc/100);
abbirb_unc.links(2).m = m2 * (1 + unc/100);
abbirb_unc.links(3).m = m3 * (1 + unc/100);
abbirb_unc.links(4).m = m456 * (1 + unc/100);
abbirb_unc.links(5).m = m456 * (1 + unc/100);
abbirb_unc.links(6).m = m456 * (1 + unc/100);

% Centers of gravity; links 3 through 6 are considered as together
abbirb_unc.links(1).r = [x_cm1 * (1 + unc/100), y_cm1 * (1 + unc/100), z_cm1 * (1 + unc/100);];
abbirb_unc.links(2).r = [x_cm2 * (1 + unc/100), y_cm2 * (1 + unc/100), z_cm2 * (1 + unc/100)];
abbirb_unc.links(3).r = [x_cm3 * (1 + unc/100), y_cm3 * (1 + unc/100), z_cm3 * (1 + unc/100)];
abbirb_unc.links(4).r = abbirb_unc.links(3).r;
abbirb_unc.links(5).r = abbirb_unc.links(3).r;
abbirb_unc.links(6).r = abbirb_unc.links(3).r;

% Inertia matrices [g * mm^2]
abbirb_unc.links(1).I = abbirb.links(1).I *(1 + unc/100);
abbirb_unc.links(2).I = abbirb.links(2).I *(1 + unc/100);
abbirb_unc.links(3).I = abbirb.links(3).I *(1 + unc/100);
abbirb_unc.links(4).I = abbirb_unc.links(3).I;
abbirb_unc.links(5).I = abbirb_unc.links(3).I;
abbirb_unc.links(6).I = abbirb_unc.links(3).I;

%% Numerical preparations

n_joint = size(abbirb.links, 2);

% Time definitions, in [s]
t_init = 0;
t_end = 10;
delta_t = 0.1;
t = t_init:delta_t:t_end;

%% Robot init

q0 = [0, 0, 0, pi/2, 0, -pi/2];
q_dot0	= [0, 0, 0, 0, 0, 0]';

q0_fig = figure2('Name', 'Initial pose');
abbirb.plot(q0);

%% Forward kinematics

% Homogeneous transformation, given the configuration
fk = abbirb.fkine(q0);

%% Trajectory

% Circumference parameters; generated based on the initial EE position
radius = 300; % [mm]
center = indexAt(fk.T, 1:3, 4) - [0; 0; radius]; % [mm]

y_traj = center(2) + radius * sin(t/t(end)*2*pi);
x_traj = center(1) * ones(size(y_traj));
z_traj = center(3) + radius * cos(t/t(end)*2*pi);

traj_fig = figure2('Name', 'Trajectory');
plot3(x_traj, y_traj, z_traj);
grid on
hold on
plot3(center(1), center(2), center(3), 'ob')
legend("trajectory", "Center (C)");
ylabel("Y [mm]");
xlabel("X [mm]");
zlabel("Z [mm]");

tr_mesh = figure2('Name', 'Trajectory and IRB1660 model with mesh');
show(abbirb_urdf, 'visuals', 'on', 'collision', 'off');
grid on
hold on
plot3(x_traj/1000, y_traj/1000, z_traj/1000);
plot3(center(1)/1000, center(2)/1000, center(3)/1000, 'ob')
xlim([-0.5, 1.2]);
ylim([-0.5, 0.5]);
zlim([0, 1.1]);
ylabel("Y [m]");
xlabel("X [m]");
zlabel("Z [m]");

sit_fig = figure2('Name', 'Initial situation');
abbirb.plot(q0, 'floorlevel', 0);
hold on
grid on
plot3(x_traj, y_traj, z_traj);
zlim([0, 1200]);
ylabel("Y [mm]");
xlabel("X [mm]");
zlabel("Z [mm]");


%% Desired configuration

% Tangent, normal and binormal to the circle, using Frenet-Sarret formulas
% (starting from the derivative of the parametric circumference used)
syms t_var real positive
tan_dir = @(t) [0, 300*2*pi/t_end*cos(t/t_end * 2*pi), -300*2*pi/t_end*sin(t/t_end*2*pi)];
tan_versor = @(t) tan_dir(t)/norm(tan_dir(t));
normal_dir = matlabFunction(diff(tan_versor(t_var)));
normal_versor = @(t) normal_dir(t)/norm(normal_dir(t));
binormal_versor = @(t) cross(tan_versor(t), normal_versor(t));

% Desired orientation based on the above versors
des_rotmat = @(t) [normal_versor(t)', tan_versor(t)', -binormal_versor(t)'];
phi_t = @(t) indexAt(MatToEulZYX(des_rotmat(t)), 1);
theta_t = @(t) indexAt(MatToEulZYX(des_rotmat(t)), 2);
psi_t = @(t) indexAt(MatToEulZYX(des_rotmat(t)), 3);

theta_des = [];
phi_des = [];
psi_des = [];
for i = t
    theta_des = [theta_des, theta_t(i)];
    phi_des = [phi_des, phi_t(i)];
    psi_des = [psi_des, psi_t(i)];
end

% Desired pose for every point of the circle
des_pose = [x_traj', y_traj', z_traj', psi_des', theta_des', phi_des'];
des_twist = [diff(des_pose,1,1); zeros(1,6)];

% Compute desired configuration via IK solution
q_des = desiredConfiguration(q0, des_pose, abbirb);
dq_des = [diff(q_des,1,2)/delta_t, zeros(6,1)];
ddq_des = [zeros(6,1), diff(dq_des,2,2)/delta_t, zeros(6,1)];

% des_dif = figure2('Name', 'Desired behavior');
% for i = 1:size(q_des, 2)
%     abbirb.plot(q_des(:, i)')
%     zlim([0, 1200]);
% end
