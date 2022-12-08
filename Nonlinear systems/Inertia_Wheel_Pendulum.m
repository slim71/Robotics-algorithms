%% Variables and constants

% NL state variables
syms theta [1 2] real
syms thetad [1 2] real

% Inertia matrices
syms J1 J2 positive real

% Inputs
syms u v

% Masses
syms mc m11 m12 m2 positive real

% Dimensions
syms l1 l2 positive real

% Friction model function and coefficients
syms fp
syms cp1 cp2 real

% System state array
state = [theta1 thetad1 theta2 thetad2];

% Parameters array
pars = [J1 J2 mc m11 m12 m2 l1 l2 cp1 cp2];

% Gravitational accelaration
g = 9.81;

% Parameters numerical values
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
first_state = [pi/4 0 0 0]; % Used in simulink

%% Vector fields and distributions describing the system

D = (m11 * l1^2 + (m12 + m2) * l2^2 + J1) * J2;
D_num = double(subs(D, pars, pars_num));
a11 = J2 / D;
a11_num = double(subs(a11, pars, pars_num));
a12 = -J2 / D;
a12_num = double(subs(a12, pars, pars_num));
a21 = a12;
a21_num = a12_num;
a22 = (m11 * l1^2 + (m12 + m2)*l2^2 + J1 + J2) / D;
a22_num = double(subs(a22, pars, pars_num));

% Control input field
g1 = [0 a12 0 a22]';
G = g1;

% Drift vector field
f = [thetad1;
    a11 * mc * g * sin(theta1) + a11;
    thetad2;
    a21 * mc * g * sin(theta1) + a21
    ];
f_wfriction = [thetad1;
    a11 * mc * g * sin(theta1) - a11 * fp;
    thetad2;
    a21 * mc * g * sin(theta1) - a21 * fp
    ];

% Distributions
delta0 = G;
delta = [f G];

% Output functions
chosen_output = a22 * theta1 + a11 * theta2;
output_wfriction = theta1;

%% Simulink environment

% Matrices defining the linearized system
A = [[0 1 0 0]; 
     [0 0 1 0]; 
     [0 0 0 1]; 
     [0 0 0 0]];
eig(A)
B = [0 0 0 1]';
% Since we've used Pole Placement to stabilize the system, here's the
% desired poles
P = [-1 -0.75 -0.5 -0.25];
Kpp = place(A, B, P);
eig(A-B*Kpp)

%% To be checked and perfectioned if needed
% Here we wanted to perform reference tracking in Simulink
% C = [0 0 0 1]; % only xi3
% A_aug = [A zeros(size(A,1),1); C 0];
% B_aug = [B; 0];
% B_ref = [zeros(rank(A),1); 1];
% P_aug = [P -10];
% K_aug = place(A_aug, B_aug, P_aug);