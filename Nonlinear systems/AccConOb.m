%% Setting up
clc; clear all; close all;

Inertia_Wheel_Pendulum

%% Small time local accessibility
% To study the stla of the system we use the Chow Filtration theorem, based on 
% the accessibility distribution

fprintf("SMALL TIME LOCAL ACCESSIBILITY\n");

delta_min = AccFiltration(delta, delta0, state);

% rank() is not accurate enough, since it doesn't account for particular 
% cases in which a symbolic matrix can lose rank. So we'll double check 
% that the determinant of the matrix is not zero when rank() returns the
% max possible rank of the considered matrix
acc_sym_rank = rank(delta_min);
acc_sym_det = det(delta_min);
fprintf("\nAs we can check from its determinant, which is \n" + ...
    "\tdet(Delta_min) = %s \nthe Accessibility matrix has \n" + ...
    "\trank(Delta_min) = %d\n", vpa(acc_sym_det), acc_sym_rank);

dmin_num = subs(delta_min, [state, pars], [state0, pars_num]); 
acc_num_rank = rank(dmin_num);
acc_num_det = det(dmin_num);

if ((acc_sym_rank == length(state)) && (double(acc_num_det) ~= 0))
    fprintf("\nThe system is stla. \n")
    is_stla = true;
else
    fprintf("\nThe system is not stla. \n")
end

%% Weak local accessibility

fprintf("\nWEAK LOCAL ACCESSIBILITY\n");

if is_stla
    fprintf("\nThe system is stla, so it is also weak locally accessible \n")
else
    delta_weak = AccFiltration(delta, delta, state);

    % see note above about rank()
    locacc_sym_rank = rank(delta_weak);
    locacc_sym_det = det(delta_weak);

    dwmin_num = subs(delta_weak, [state, pars], [state0, pars_num]);
    locacc_num_rank = rank(dwmin_num);
    locacc_num_det = det(dwmin_num);

    if ((locacc_sym_rank == length(state)) && (doouble(locacc_num_det) ~= 0))
        fprintf("\nThe system is wla. \n")
        is_stla = true;
    else
        fprintf("\nThe system is not wla. \n")
    end
    
end

%% Small time local controllability

fprintf("\nSMALL TIME LOCAL CONTROLLABILITY\n");

[is_stlc, matched_cond] = STLCCheck(f, G, state, pars, state0, pars_num);

if is_stlc
    fprintf("\nThe Local Controllability condition %d was met, so the " + ...
            "system is stlc. \n", matched_cond);
    is_stlc = true;
else
    fprintf("\nThe system is not stlc \n")
end
%% Local observability

fprintf("\nLOCAL OBSERVABILITY\n");

h = chosen_output;
fprintf("\nAs output of the system we choose: \n h = %s\n", h)

dh = jacobian(h, state);
omega0 = dh;

omega_max = ObsFiltration(delta, omega0, state);

% see note above about rank()
obs_sym_rank = rank(omega_max);
obs_sym_det = det(omega_max);
fprintf("\nAs we can check from its determinant, which is \n" + ...
    "\tdet(\xOmega_max) = %s \nthe Observability matrix has \n" + ...
    "\trank(Omega_max) = %d\n", vpa(obs_sym_det), obs_sym_rank);

omax_num = subs(omega_max, [state, pars], [state0, pars_num]); 
obs_num_rank = rank(omax_num);
obs_num_det = det(omax_num);

if ((obs_sym_rank == length(state)) && (double(obs_num_det) ~= 0))
    fprintf("\nThe Observability matrix has full rank in the general symbolic " + ...
        "form, so the system can be locally observable with this output " + ...
        "choice. \n");

    if obs_num_rank == length(state) % we already checked the determinant
        fprintf("\nThe Observability matrix has full rank in the specified " + ...
            "inital conditions, so the system is locally observable \n");
    else
        fprintf("\nThe Observability matrix has rank %d in the specified inital " + ...
            "conditions, so the system is not locally observable \n", obs_num_rank);
    end

else
    fprintf("\nThe Observability matrix has rank %d in the general symbolic " + ...
        "form, so the system cannot be locally observable with this output " + ...
        "choice. \n", obs_sym_rank);
end
