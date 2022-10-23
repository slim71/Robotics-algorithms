%% Setting up
clc; clear all; close all;
digits 4;

Inertia_Wheel_Pendulum

F_nopars = subs(f, pars, pars_num);
F_nopars_fric = subs(f_wfriction, pars, pars_num);
G_nopars = subs(G, pars, pars_num);

%% Complete feedback linearization
fprintf("\nCOMPLETE FEEDBACK LINEARIZATION\n");

[rel_deg, Lf_total, Lg_total] = SISOFL(f, G, chosen_output, state, state0);

fprintf("\nThe system has relative degree %d and the state space has " + ...
    "size %d.\n", rel_deg, length(state));

if rel_deg >= length(state)
    fprintf("\nSo we can achieve a complete decoupling feedback " + ...
        "linearization. \n")
else
    fprintf("\nSo a complete decoupling feedback linearization is not " + ...
        "possible. We must proceed in another way. \n")
end

xi_value = Lf_total(1:end-1);
alpha = - Lf_total(end) / Lg_total(end);
beta = 1 / Lg_total(end);

fprintf("\nWe can perform the change of variables \n \x3A6(x) = \n");
fprintf("\t\t%s\n", xi_value);

u = alpha + beta * v;
fprintf("\nThe feedback linearization control we look for is then \nu = %s\n", u);
fprintf("and following the above results we have y^(r) = v \n");

syms xi [rel_deg 1]
syms xi_dot [rel_deg 1]
syms y_dotr

new_f = [xi_value(2:end).'; Lf_total(end)];
new_g = [zeros(rel_deg-1, 1); Lg_total(end)];
v = y_dotr;

% xi_dot = new_f + new_g * v
xi_dot(1) = new_f(1) + new_g(1) * v;
xi_dot(2) = new_f(2) + new_g(2) * v;
xi_dot(3) = new_f(3) + new_g(3) * v;
xi_dot(4) = new_f(4) + new_g(4) * v;
fprintf("\nWe then have the new system: \n \x3BE = \n");
fprintf("\t%s\n", xi_dot);

fprintf("\nFrom here on out, we could add another feedback control to " + ...
    "stabilize the system and give it the eigenvalues we want, \nor we " + ...
    "could also design a particular controller for what we're trying " + ...
    "to achieve. \n");

%% Partial feedback linearization
fprintf("\nPARTIAL FEEDBACK LINEARIZATION\n");

[r_fric, Lf_fric, Lg_fric] = SISOFL(f, G, output_wfriction, state, state0);

fprintf("\nConsidering frictions in the output, the system has " + ...
    "relative degree %d; the state space is still of " + ...
    "size %d.\n", r_fric, length(state));

if r_fric >= length(state)
    fprintf("\nWe can still achieve a complete decoupling feedback " + ...
        "linearization. \n")
else
    fprintf("\nNow a complete decoupling feedback linearization is not " + ...
        "possible. We must proceed in another way. \n")
end

xi_fric = Lf_fric(1:end-1);
zeta1 = a11 * fp + theta2;
zeta2 = thetad2;
zeta = [zeta1, zeta2];
fprintf("\nWe first choose as state variables \n\x3BE = \n");
fprintf("\t%s\n", xi_fric);
fprintf("\nan to complete the variable change we set: \n \x3B6 = \n");
fprintf("\t%s\n", zeta);

fprintf("\nIt can be easily checked that these are indipendent from " + ...
    "the \x3BE variables already established, and among themselves. \n" + ...
    "We have \n\trank([[\x3BE], [\x3B6]])=%d \nand \n\t" + ...
    "det(\x3A6)|_0=det([\x3BE; \x3B6])|_0=%d \n", ...
    rank([xi_fric.',[theta1+theta2; theta2]]), ...
    subs(det(jacobian([xi_fric, zeta], state)), state, state0));

Lg_zeta1 = ScalarFDerivative(G, zeta1, state);
Lg_zeta2 = ScalarFDerivative(G, zeta2, state);
fprintf("\nMoreover: \n Lg\x3B6\x2081 = %s \n Lg\x3B6\x2082 = %s \n", ...
    Lg_zeta1, Lg_zeta2)

syms dxi_pfl [1 2]
xi_pfl = [xi(1); xi(2)];
dxi_pfl(1) = Lf_fric(2);
dxi_pfl(2) = Lf_fric(end) + Lg_fric(end) * u;
fprintf("\nAfter the change of variables \x3A6(x) = [\x3BE; \x3B6], " + ...
    "we have \nd\x3BE = \n");
fprintf("\t%s\n", dxi_pfl);

a = Lg_fric(end);
b = Lf_fric(end);
u = - b/a + 1/a * v;
fprintf("\nFor the partial feedback linearization we then set: v = \x3BE^(r)\n");
fprintf("as new input and the feedback control as: \n u = %s\n", u);

% Compute f-derivatives in x_0
Lf_zeta1 = ScalarFDerivative(f, zeta1, state);
Lf_zeta2 = ScalarFDerivative(f, zeta2, state);
Lf_zeta1_0 = subs(Lf_zeta1, xi_fric, [0, 0]);
Lf_zeta2_0 = subs(Lf_zeta2, xi_fric, [0, 0]);

% Linearized system dynamics
r = r_fric; % = length(xi_fric)
n = length(state);
m = length(zeta);
q = [Lf_zeta1; Lf_zeta2];
p = [Lg_zeta1; Lg_zeta2];
A = [[zeros(r-1, 1), eye(r-1)]; zeros(1, r)];
B = [zeros(r-1, 1); 1];
C = [1 zeros(1, r-1)];

fprintf("\nThe linearized system would be: \n");
dynamics = A*xi_pfl + B*v
out = C*xi_pfl

%% 0-dynamics
fprintf("\n ZERO DYNAMICS\n");

fprintf("\nAfter achieving a partial feedback linearization, we want " + ...
    "to make sure the linearized system is also internally stable. \n");

fprintf("\nTo ensure y(t)=0, d\x3B6 = q(0, \x3B6) should be " + ...
    "asympthotically stable in \x3B6_0.\n");

q0 = subs(q, xi_fric, [0,0]);
p0 = subs(p, xi_fric, [0,0]);
zeta1_dot = q0(1) + p0(1)*v
zeta2_dot = q0(2) + p0(2)*v

q_mat = [[0, 1]; [-1/(fp + theta2*(m11*l1^2 + (m2 + m12)*l2^2 + J1)), 0]];
q_svd = svd(q_mat);
fprintf("\nNow, q(\x3BE\x2080) can be written as\nq(\x3BE\x2080) = \n");
fprintf("\t\t%s %s\n", q_mat.');
fprintf("which has singular values \n\x3BB = \n");
fprintf("\t%s\n", q_svd);

friction = @(th1d) cp1*sign(th1d) + cp2*th1d;
fprintf("\nFinally, the friction can be modelled as %s\n", char(friction));

svd_num = subs(q_svd, pars, pars_num);
fprintf("\nSince d\x3B8\x2081 is \x3Be\x2082, for the zero dynamics " + ...
    "we must set it to zero, obtaining \n\x3BB = \n")
fprintf("\t%s\n", vpa(subs(svd_num, fp, friction(0))));

fprintf("\nSince these are always positive for any 𝜁_0, the zero " + ...
    "dynamics is not stable and we cannot conclude on the nonlinear " + ...
    "system internal stability.\n");