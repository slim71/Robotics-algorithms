function qsol = sot(J, xdot, kind, lambda)
% SOT resolves the Stack Of tasks problem and return the join coordinates.
%
% Input
%   J:      vector containing the tasks Jacobians
%   xdot:   vector containing the 
% Output
%   qdot:   
%
    syms t
    syms q [6 1]
    
    time = 2*pi;
    step = 0.1;
    n_task = size(xdot,2) + 1; % == size(J,3)
    n_joints = size(J, 2);
%     qdot = zeros(n_joints, n_task);
    qdot = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01]/1000;
    PA = eye(n_joints);

    for j = [1:n_task-1]
        PA = cat(3, PA, eye(n_joints));
    end
    
    qapp = [];
    
    % Extend xdot and J to work with starting index 1
    % to account for the previous instant in the algorithm
    J = cat(3, zeros(size(J,1),size(J,2)), J);
    xdot = [zeros(size(xdot,1), 1), xdot];
    
    for timestep = [0:step:time]
        tt = timestep;
        
        if (kind == 0)
            for i = [2:n_task]
                qdot(:, i) = qdot(:, i-1) + pinv(subs(J(:, :, i), q, qdot(:, i-1)) * PA(:, :, i-1)) * (subs(xdot(:, i), t, tt) - subs(J(:, :, i), q, qdot(:, i-1)) * qdot(:, i-1));
                PA(:, :, i) = PA(:, :, i-1) - pinv(subs(J(:, :, i), q, qdot(:, i-1)) * PA(:,:,i-1)) * subs(J(:, :, i), q, qdot(:, i-1)) * PA(:,:,i-1);
            end
        elseif (kind == 1)
            for i = [2:n_task]
                qdot(:, i) = qdot(:, i-1) + damped(J(:, :, i) * PA(:, :, i-1), lambda) * (xdot(:, i) - J(:, :, i) * qdot(:, i-1));
                PA(:, :, i) = PA(:, :, i-1) - damped(J(:,:,i) * PA(:,:,i-1), lambda) * J(:,:,i) * PA(:,:,i-1);
            end
        end
        
        qapp = cat(3, qapp, qdot);
    end
    
    qsol = qapp;%qdot(:,2:end);
end