% [U,S,V] = svd(A);
% n: dof, joint dim
% m: task dim
% m<=n, sum mk<=n
% J: mxn, full rank rho=m
% #: pseudo Moore-Penrose
% P=I-J^#*J: nxn
% qndot0=0
% i=1...k-1
% k=1...l

% TODO: function helps
% TODO: check type and dimension of arguments

affirmation = {'not', ''};

%% 1st example from the draft
fprintf("\t First example \n")

J1 = [1 0 0; 0 1 0];
J2 = [1 0 0];

fprintf("J1 and J2 are %s independent. \n\n", ...
    affirmation{check_independence(J1, J2)+1})

fprintf("Obtained first jacobian projector: \n ")
P1_1 = null_projector(J1);
fprintf("\t%f %f %f \n", P1_1.')

[ind1, num_cond1] = check_ill_conditioning(J2*P1_1);
fprintf("\nJ2*P1 is %s ill-conditioned. \n\t J2*P1 \t\t pseudoinverse: \n",...
    affirmation{ind1+1})
pseudo1 = pinv(J2*P1_1);
fprintf("\t%f\t%f \n", J2*P1_1, pseudo1)

%% Second example
fprintf("\n\t Second example \n")

J1 = [1 0 0; 0 1 0];
J2 = [1 0 0.000001];

fprintf("J1 and J2 are %s independent. \n\n", ...
    affirmation{check_independence(J1, J2)+1})

fprintf("Obtained first jacobian projector: \n")
P1_2 = null_projector(J1);
fprintf("\t%f %f %f \n", P1_2.')

[ind2, num_cond2] = check_ill_conditioning(J2*P1_2);
fprintf("\nJ2*P1 is %s ill-conditioned. \n \t J2*P1 \t\t pseudoinverse: \n",...
    affirmation{ind2+1})
pseudo2 = pinv(J2*P1_2);
fprintf("\t%f\t%f \n", J2*P1_2, pseudo2)


% cond_num = norm(pinv(matrix))*norm(matrix)

%% Trying damped pseudoinverse

A =[[0.5705,    0     -0.8213        0     ]
     [0.8213,    0      0.5705        0     ]
     [0,        -1.0    0           -11.0000]
     [0,         0      0             1.0   ]];
 
 invA = inv(A);
 sigma_invA = svd(invA);
 
 pA = pinv(A);
 sigma_pA = svd(pA);
 
 dampedA = damped(A, 100);
 sigma_dA = svd(dampedA);