%% Setting up
clc; clear all; close all;

Inertia_Wheel_Pendulum

%% Small time local accessibility
% To study the stla of the system we use the Chow Filtration theorem, based on 
% the accessibility distribution

fprintf("SMALL TIME LOCAL ACCESSIBILITY\n");

delta_min = AccFiltration(delta, delta0, state);

if rank(subs(delta_min, [state, pars], [state0, pars_num])) == length(state)
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

    if rank(subs(delta_weak, [state, pars], [state0, pars_num])) == length(state)
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

rgm = rank(omega_max);

if rgm == length(state)
    fprintf("\nThe Observability matrix has full rank in the general symbolic " + ...
        "form, so the system can be locally observable with this output " + ...
        "choice. \n");

    gamma_min_num = subs(omega_max, [state, pars], [state0, pars_num]);
    rgm_num = rank(gamma_min_num);

    if rgm_num == length(state)
        fprintf("\nThe Observability matrix has full rank in the specified " + ...
            "inital conditions, so the system is locally observable \n");
    else
        fprintf("\nThe Observability matrix has rank %d in the specified inital " + ...
            "conditions, so the system is not locally observable \n", rgm_num);
    end

else
    fprintf("\nThe Observability matrix has rank %d in the general symbolic " + ...
        "form, so the system cannot be locally observable with this output " + ...
        "choice. \n", rgm);
end
