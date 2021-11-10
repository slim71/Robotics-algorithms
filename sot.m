function qsol = sot(J, xdot, kind, lambda)
% SOT resolves the Stack Of Tasks problem and return the join coordinates.
%
%     xdot: 6 x n_tasks
%     J:    (6 x n_joints) x n_tasks
%
%     qdot: (n_joints x n_tasks) x n_steps
%     PA:   (n_joints x 6) x n_tasks

    syms t
    % TODO: make time variables as arguments?
    total_time = 2*pi;
    step = 0.01;
    n_steps = ceil(total_time/step);
    n_joints = size(J, 2);
    n_tasks = size(xdot,2);
    syms q [n_joints 1] % TODO: automatically get sym variable?

    % +1 because the 1-index column of all pages will be qdot0(t)=0
    qdot = zeros(n_joints, n_tasks+1, n_steps); %checked: OK to be 0s
%     qdot(:,:,1) = repmat(pi/10, n_joints, n_tasks+1);
    % Resulting in a (n_joints*1 x n_joints*1 x 1*n_tasks) matrix
    PA = repmat( eye(n_joints), 1, 1, n_tasks+1); % checked: OK to be an Identity
    % Extends twist vector to use the same index as the other matrices
    twists = cat(2,zeros(size(xdot,1),1), xdot);

    js = cat(3, zeros(size(J,1),size(J,2)),J);


%     for timestep = [step:step:total_time] % or [step:step:n_steps*step]
    for i=[1:n_tasks] % i=2 <-> task 1
        current_task = i+1;
        prev_task = i;
        
        proj_prev = double(PA(:, :, prev_task)); % Projector of J_i-1

        for curr_time_step = [1:n_steps]
            timestep = curr_time_step * step;
            if curr_time_step <= 1
                prev_step = curr_time_step;
            else
                prev_step = curr_time_step-1;
            end

            % Jacobian: current task, previous timestep approximation
            jac = double(subs(js(:,:,current_task), q, qdot(:, current_task, prev_step)));

            % Desired twist: current task, current timestep
            xx = double(subs(twists(:, current_task), t, timestep));

            % Previous task final solution
            q_prev = qdot(:, prev_task, end);
    
            % Compute next timestep solution, for each task i
            q_next = q_prev + pinv(jac * proj_prev) * (xx - jac * q_prev);

            % i-th task
            qdot(:, current_task, curr_time_step) = q_next;
        end

        proj_next = proj_prev - pinv(jac * proj_prev) * jac * proj_prev;
        PA(:, :, current_task) = proj_next;
    end


    qsol = qdot(:,2:end,:);
end