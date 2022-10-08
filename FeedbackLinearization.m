%% Setting up
clc; clear all; close all;
digits 4;

Inertia_Wheel_Pendulum

F_nopars = subs(f, pars, pars_num);
F_nopars_fric = subs(f_wfriction, pars, pars_num);
G_nopars = subs(G, pars, pars_num);

%% I/O feedback linearization analysis

[rel_deg, Lf_total, Lg_total] = SISOFL(F_nopars, G_nopars, ...
                       subs(chosen_output, pars, pars_num), state, state0);

fprintf("The system has relative degree %d and the state space has " + ...
    "size %d.\n", rel_deg, length(state));

if rel_deg >= length(state)
    fprintf("So we can achieve a complete decoupling feedback " + ...
        "linearization. \n")
else
    fprintf("So a complete decoupling feedback linearization is not " + ...
        "possible. We must proceed in another way. \n")
end

%% Complete feedback linearization

xi_value = Lf_total(1:end-1);
alpha = - Lf_total(end) / Lg_total(end);
beta = 1 / Lg_total(end);

fprintf("We can perform the change of variables phi(x) = ");
xi_value.'

fprintf("The feedback linearization we look for is then");
u_fl = vpa(alpha + beta * v)
fprintf("and following the above results we have y^(r) = v \n");

fprintf("We then have the new system: \n");
syms xi [rel_deg 1]
syms xi_dot [rel_deg 1]
syms doty_r
new_state = xi;
new_state_dot = xi_dot;
new_f = [vpa(xi_value(2:end)).'; vpa(Lf_total(end))];
new_g = [zeros(rel_deg-1, 1); vpa(Lg_total(end))];
% xi_dot = new_f + new_g * v
xi_dot1 = new_f(1) + new_g(1) * doty_r
xi_dot2 = new_f(2) + new_g(2) * doty_r
xi_dot3 = new_f(3) + new_g(3) * doty_r
xi_dot4 = new_f(4) + new_g(4) * doty_r

fprintf("From here on out, we could add another feedback control to " + ...
    "stabilize the system and give it the eigenvalues we want, \nor we " + ...
    "could also design a particular controller for what we're trying " + ...
    "to achieve. \n");

%% Partial feedback linearization

fprintf("\n");

[r_fric, Lf_fric, Lg_fric] = SISOFL(F_nopars_fric, G_nopars, ...
                       subs(output_wfriction, pars, pars_num), state, state0);

fprintf("Considering frictions in the output, the system has " + ...
    "relative degree %d; the state space is still of " + ...
    "size %d.\n", r_fric, length(state));

if r_fric >= length(state)
    fprintf("We can still achieve a complete decoupling feedback " + ...
        "linearization. \n")
else
    fprintf("Now a complete decoupling feedback linearization is not " + ...
        "possible. We must proceed in another way. \n")
end

xi_fric = Lf_fric(1:end-1);

fprintf("To complete the variable change we choose: \n");

syms zeta1 zeta2
zeta1_fric = a11 * fp + theta2;
zeta2_fric = thetad2;
zeta_fric = [zeta1_fric, zeta2_fric];

fprintf("It can be easily checked that these are indipendent from " + ...
    "the xi variables already established, and among themselves. \n" + ...
    "We have \nrank([[xi1;x2], [zeta1, zeta2]])=%d \n and \n" + ...
    "det(phi)|_0=det([xi; zeta])|_0=%d \n", ...
    rank([xi_fric.',[theta1+theta2; theta2]]), ...
    subs(det(jacobian([xi_fric, zeta_fric], state)), state, state0));

fprintf("Moreover: \n")
Lg_zeta1 = ScalarFDerivative(G_nopars, zeta1_fric, state)
Lg_zeta2 = ScalarFDerivative(G_nopars, zeta2_fric, state)

phi = [xi_fric, zeta_fric];

fprintf("After the change of variables wee have \n");
syms xi_dot [1 2]
xis = [xi(1); xi(2)];
xi1_dot = Lf_fric(2)
xi2_dot = Lf_fric(end) + Lg_fric(end) * u

a = Lg_fric(end);
b = Lf_fric(end);

fprintf("For the partial feedback linearization we then set: \n");
xi_dot2 = v
fprintf("as new input and the feedback control as: \n");
u = - b/a + 1/a * v

Lf_zeta1 = ScalarFDerivative(F_nopars_fric, zeta1_fric, state);
Lf_zeta2 = ScalarFDerivative(F_nopars_fric, zeta2_fric, state);
Lf_zeta1_0 = subs(Lf_zeta1, xi_fric, [0, 0]);
Lf_zeta2_0 = subs(Lf_zeta2, xi_fric, [0, 0]);

r = r_fric; % = length(xi_fric)
n = length(state);
m = length(zeta_fric);
q = [Lf_zeta1; Lf_zeta2];
p = [Lg_zeta1; Lg_zeta2];
A = [[zeros(r-1, 1), eye(r-1)]; zeros(1, r)];
B = [zeros(r-1, 1); 1];
C = [1 zeros(1, r-1)];

% Linearized system
dyn = A*xis + B*v;
out = C*xis;

%% 0-dynamics
fprintf("After achieving a partial feedback linearization, we must " + ...
    "make sure the linearized system is also internally stable. \n");

fprintf("To ensure y(t)=0 it must be zeta_dot = q(0, eta) for all zeta_0\n");

q0 = subs(q, xi_fric, [0,0]);
zeta1_dot = q0(1) + p(1)
zeta2_dot = q0(2) + p(2)

fprintf("Since this is not the case, the feedback linearization " + ...
    "control does not make the non-linear system both internally " + ...
    "and externally stable. \nNo conclusion can be made with only this " + ...
    "theorem.\n");