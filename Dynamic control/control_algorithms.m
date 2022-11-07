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

% robot_fig = figure2('Name', 'Robot model');
% show(abbirb_urdf);
% xlim([-0.5, 1]);
% ylim([-0.5, 0.5]);
% zlim([0, 1.1]);

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
abbirb.links(4).m = 0;
abbirb.links(5).m = 0;
abbirb.links(6).m = 0;

% Centers of gravity; links 3 through 6 are considered as together
abbirb.links(1).r = [x_cm1, y_cm1, z_cm1];
abbirb.links(2).r = [x_cm2, y_cm2, z_cm2];
abbirb.links(3).r = [x_cm3, y_cm3, z_cm3];
abbirb.links(4).r = [0, 0, 0];
abbirb.links(5).r = [0, 0, 0];
abbirb.links(6).r = [0, 0, 0];

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
abbirb.links(4).I = diag([0, 0, 0]);
abbirb.links(5).I = diag([0, 0, 0]);
abbirb.links(6).I = diag([0, 0, 0]);

%% Uncertain model parameters

% Percentage of uncertainty on parameters value
unc = 10; % 10%

% Masses
abbirb_unc.links(1).m = m1 * (1 + unc/100);
abbirb_unc.links(2).m = m2 * (1 + unc/100);
abbirb_unc.links(3).m = m3 * (1 + unc/100);
abbirb_unc.links(4).m = 0;
abbirb_unc.links(5).m = 0;
abbirb_unc.links(6).m = 0;

% Centers of gravity; links 3 through 6 are considered as together
abbirb_unc.links(1).r = [x_cm1 * (1 + unc/100), y_cm1 * (1 + unc/100), z_cm1 * (1 + unc/100);];
abbirb_unc.links(2).r = [x_cm2 * (1 + unc/100), y_cm2 * (1 + unc/100), z_cm2 * (1 + unc/100)];
abbirb_unc.links(3).r = [x_cm3 * (1 + unc/100), y_cm3 * (1 + unc/100), z_cm3 * (1 + unc/100)];
abbirb_unc.links(4).r = [0, 0, 0];
abbirb_unc.links(5).r = [0, 0, 0];
abbirb_unc.links(6).r = [0, 0, 0];

% Inertia matrices [g * mm^2]
abbirb_unc.links(1).I = abbirb.links(1).I *(1 + unc/100);
abbirb_unc.links(2).I = abbirb.links(2).I *(1 + unc/100);
abbirb_unc.links(3).I = abbirb.links(3).I *(1 + unc/100);
abbirb_unc.links(4).I = diag([0, 0, 0]);
abbirb_unc.links(5).I = diag([0, 0, 0]);
abbirb_unc.links(6).I = diag([0, 0, 0]);

%% Numerical preparations

n_joint = size(abbirb.links, 2);

% Time definitions, in [s]
% TODO: check, maybe to change
t_init = 0;
t_end = 10;
delta_t = 0.1;
t = t_init:delta_t:t_end;

%% Robot init

q0 = [0, 0, 0, 0, 0, 0];
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

x_traj = center(1) + radius * sin(t/t(end)*2*pi);
y_traj = center(2) * ones(size(x_traj));
z_traj = center(3) + radius * cos(t/t(end)*2*pi);

traj_fig = figure2('Name', 'Trajectory');
plot3(x_traj, y_traj, z_traj);
grid on
hold on
plot3(center(1), center(2), center(3), 'ob')
legend("trajectory", "Center (C)")
ylabel("Y [mm]");
xlabel("X [mm]");
zlabel("Z [mm]");

sit_fig = figure2('Name', 'Initial situation');
hold on
plot3(x_traj, y_traj, z_traj);
abbirb.plot(q0);
grid on

%% Desired configuration

angles = MatToEulZYX(indexAt(fk.T, 1:3, 1:3));
pose_size = size(x_traj, 2);

% Desired orientation based on the initial one
des_theta = repmat(angles(1), 1, pose_size);
des_phi = repmat(angles(2), 1, pose_size);
des_psi	= repmat(angles(3), 1, pose_size);

% Desired pose for every point of the circle
des_pose = [x_traj', y_traj', z_traj', des_theta', des_phi' des_psi'];

% Compute desired configuration via IK solution
q_des = desiredConfiguration(q0, des_pose, abbirb);
dq_des	= gradient(q_des)/delta_t;
ddq_des = gradient(dq_des)/delta_t;

des_dif = figure2('Name', 'Desired behavior');
for i = 1:size(q_des, 2)
    abbirb.plot(q_des(:, i)')
end
