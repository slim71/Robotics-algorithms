%% Variables and constants

syms theta [1 2] real

syms thetad [1 2] real

% Inertia matrix
syms J1 J2 positive real

% tau
syms u

%  mass
syms mc m11 m12 m2 positive real

% 
syms l [1 2] positive real

syms cp1 cp2 real

% System state array
state = [theta1 thetad1 theta2 thetad2];

% Parameters array
pars = [J1 J2 mc m11 m12 m2 l1 l2 cp1 cp2];

% Gravitational accelaration
g = 9.81;

% Numerical values out of "Robotics, vision and control - P. Corke" example
J1_num = 180E-6; % [kg/m^2]
J2_num = 36E-6; % [kg/m^2]
mc_num = 4.34; % [kg]
m11_num = 0.03; % [kg]
m12_num = 0.06; % [kg]
m2_num = 0.235; % [kg]
l1_num = 0.04; % [m]
l2_num = 0.125; % [m]
cp1_num = 0.00225; % [(N*m*s)/rad]
cp2_num = 0.0002; % [(N*m*s)/rad]

pars_num = [J1_num J2_num mc_num m11_num m12_num m2_num l1_num l2_num cp1_num cp2_num];

% Flags
is_stla = false;
is_stlc = false;
is_lo = false;

%% Initial conditions
state0 = [0 0 0 0];

%% Dynamic system
fprintf("System considered:")
% Xdot = w * (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) - v * (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta)) + u * (cos(psi) * cos(theta))
% Ydot = v * (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) - w * (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta)) + u * (cos(theta) * sin(psi))
% Zdot = w * (cos(phi) * cos(theta)) - u * (sin(theta)) + v * (cos(theta) * sin(phi))
% Phidot = p + r * (cos(phi) * tan(theta)) + q * (sin(phi) * tan(theta))
% Thetadot = q * cos(phi) - r * sin(phi)
% Psidot = r * cos(phi) / cos(theta) + q * sin(phi) / cos(theta)
% Pdot = (Iy - Iz) / Ix * r * q + (tau_x + tau_wx) / Ix
% Qdot = (Iz - Ix) / Iy * p * r + (tau_y + tau_wy) / Iy
% Rdot = (Ix - Iy) / Iz * p * q + (tau_z + tau_wz) / Iz
% Udot = r * v - q * w - g * sin(theta) + fw_x / m
% Vdot = p * w - r * u + g * sin(phi) * cos(theta) + fw_y / m
% Wdot = q * u - p * v + g * cos(theta) * cos(phi) + (fw_z - f_t) / m

%% Vector fields and distributions describing the system

D = (m11 * l1^2 + (m12 + m2) * l2^2 + J1) * J2;
a11 = J2 / D;
a12 = -J2 / D;
a21 = a12;
a22 = (m11 * l1^2 + (m12 + m2)*l2^2 + J1 + J2) / D;

% Control input field
g1 = [0 a12 0 a22]';
G = g1;

% Drift vector field
f = [thetad1;
    a11 * mc * g * sin(theta1) + a11;
    thetad2;
    a21 * mc * g * sin(theta1) + a21
    ];

% Distributions
delta0 = G;
delta = [f G];

% Output functions
chosen_output = a22 * theta1 + a11 * theta2;