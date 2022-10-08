%% Setting up
clc; clear all; close all;

Inertia_Wheel_Pendulum

%% Small time local accessibility
% To study the stla of the system we use the Chow Filtration theorem, based on 
% the accessibility distribution

delta_min = AccFiltration(delta, delta0, state);

if rank(subs(delta_min, [state, pars], [state0, pars_num])) == length(state)
    fprintf("The system is stla \n")
    is_stla = true;
else
    fprintf("The system is not stla \n")
end

%% Weak local accessibility

if is_stla
    fprintf("The system is stla, so it is also weak locally accessible \n")
else
    delta_weak = AccFiltration(delta, delta, state);

    if rank(subs(delta_weak, [state, pars], [state0, pars_num])) == length(state)
        fprintf("The system is wla. \n")
        is_stla = true;
    else
        fprintf("The system is not wla. \n")
    end
    
end

%% Small time local controllability

[is_stlc, matched_cond] = STLCCheck(f, G, state, pars, state0, pars_num);

if is_stlc
    fprintf("The Local Controllability condition %d was met, so the " + ...
            "system is stlc. \n", matched_cond);
    is_stlc = true;
else
    fprintf("The system is not stlc \n")
end
%% Local observability

fprintf("As output of the system we choose: \n")
h = chosen_output 

dh = jacobian(h, state);
gamma0 = dh;

gamma_min = ObsFiltration(delta, gamma0, state);

rgm = rank(gamma_min);

if rgm == length(state)
    fprintf("The Observability matrix has full rank in the general symbolic " + ...
        "form, so the system can be locally observable with this output " + ...
        "choice. \n");

    gamma_min_num = subs(gamma_min, [state, pars], [state0, pars_num]);
    rgm_num = rank(gamma_min_num);

    if rgm_num == length(state)
        fprintf("The Observability matrix has full rank in the specified " + ...
            "inital conditions, so the system is locally observable \n");
    else
        fprintf("The Observability matrix has rank %d in the specified inital " + ...
            "conditions, so the system is not locally observable \n", rgm_num);
    end

else
    fprintf("The Observability matrix has rank %d in the general symbolic " + ...
        "form, so the system cannot be locally observable with this output " + ...
        "choice. \n", rgm);
end
