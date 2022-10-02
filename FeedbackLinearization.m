%% Setting up
clc; clear all; close all;

% Quadrotor_system
Simplified_quadrotor_system

%% MIMO conditions

F_nopars = subs(f, pars, pars_num);
G_nopars = subs(G, pars, pars_num);

%% I/O feedback linearization

[r_mat, poss_epsi, Gamma_mat, E_mat] = MIMOFL(F_nopars, G_nopars, chosen_output, state, state0);

% Compute relative degree for each output
rel_degrees = zeros(size(r_mat, 1), 1);
for i = 1:size(r_mat, 1)
    r_app = r_mat(i, :);

    rd = min(r_app(r_app>0));
    if isempty(rd) % No definite relative degree for this output
        rel_degrees(i, 1) = abs(max(r_app(r_app<0)));
    else % At least a definite relative degree for this output
        rel_degrees(i, 1) = rd;
    end
end

E = subs(E_mat, state, state0);
if rank(E) == min(size(E))
    fprintf("E is full-rank, so it's invertible. \n");

    % Custom formatter for fprintf
    outputstr = repmat('%i ', 1, size(rel_degrees, 1)-1);
    outputstr = [outputstr '%i'];

    fprintf("Then its relative degree is [" + outputstr + ...
        "] (total: %d). \n", rel_degrees.', sum(rel_degrees));
else
    fprintf("E is not invertible (rank: %d, dimensions: [%d %d]). \n", ...
            rank(E), size(E));
    fprintf("So it does not have a definite relative degree.\n");
end

%% Dynamic feedback linearization
disp(E_mat)
fprintf("Looking once more at the symbolic versione of E, we see " + ...
        "that the 2nd column is comprised of all 0s.\n" + ...
        "This means the second input is in some way 'slower' than " + ...
        "the others. \n" + ...
        "A possible solution, explored here, is to slow down the " + ...
        "other inputs with one or more integrators.\n");

fprintf("Here we set the first input equal to the output of a " + ...
        "double integrator driven by the u1_bar:" + ...
        "u1 = mu;  mu_dot = nu; nu_dot = u1_bar \n");
fprintf("Now u1 is not an input anymore, but rather a part of the " + ...
        "internal state. In its place, u1_bar becomes the new input.\n");

syms zeta xi
aug_state = [state zeta xi];
aug_state0 = [state0 0 0];
aug_f = [F_nopars(1:6); g7_1; g8_1; g9_1; F_nopars(10:end); xi; 0];
g1_bar = [zeros(size(g1, 1), 1); 0; 1];
aug_g2 = subs([g2; 0; 0], pars, pars_num);
aug_g3 = subs([g3; 0; 0], pars, pars_num);
aug_g4 = subs([g4; 0; 0], pars, pars_num);
aug_G = [g1_bar aug_g2 aug_g3 aug_g4];

[aug_r, aug_poss_epsi, aug_Gamma, aug_E] = MIMOFL(aug_f, aug_G, chosen_output, aug_state, aug_state0);

% Compute new relative degree for each output
aug_rel_degrees = zeros(size(aug_r, 1), 1);
for i = 1:size(aug_r, 1)
    r_app = aug_r(i, :);

    rd = min(r_app(r_app>0));
    if isempty(rd) % No definite relative degree for this output
        aug_rel_degrees(i, 1) = abs(max(r_app(r_app<0)));
    else % At least a definite relative degree for this output
        aug_rel_degrees(i, 1) = rd;
    end
end

% aug_E is not full-rank only because of the limit imposed on the symbolic
% calculations in SISOFL [line 12] (withouth them my PC would not finish
% computing higher order derivatives with symbolic calculations

% Custom formatter for fprintf
outputstr = repmat('%i ', 1, size(aug_rel_degrees, 1)-1);
outputstr = [outputstr '%i'];

fprintf("The relative degree of the newly computed E is [" ...
        + outputstr + "] (total: %d). \n", aug_rel_degrees.', ...
        sum(aug_rel_degrees));

%%

phi = [x(1);
    x(3)-x(2)^3;
    x(2)];

det(subs(jacobian(phi, [x(1),x(2),x(3)]), [x(1),x(2),x(3)], [0,0,0]))
