function qdot = sot(tt, qq, sym_Js, sym_twistes, lambda)
% TODO help

    % Initialization
    q0dot = zeros(max(size(qq)),1);

    qdoti = q0dot;
    Pi = eye(max(size(qq)));

    for i = 1:max(size(sym_twistes))
        Ji = sym_Js{i}; % Jacobian function handle
        twisti = sym_twistes{i}; % Task function handle

        % Computed Jacobian
        Ji_qq = Ji(qq);

        [log, ~] = check_ill_conditioning(Ji_qq * Pi);
        if log
            % Using damped pseudo-inversion
            inverted = damped_pinv(Ji_qq*Pi, lambda);
        else
            % Standard pseudo inversion
            inverted = pinv(Ji_qq*Pi);
        end

        
        % Solution for i=2...n:
        % q'_i = q'_(i-1) + inv(J_i*P_(i-1)) * (task_i - J_i*q'_(i-1))
        qdoti = qdoti + inverted * (twisti(tt,qq)' - Ji_qq * qdoti);

        % Standard way to compute projector P:
        % Pi = eye(max(size(qq))) - pinv(tot_J)*tot_J; 

        % Iterative method used here:
        % P_i = P_(i-1) - inv(J_i * P_(i-1)) * J_i * P_(i-1)
        Pi = Pi - inverted * Ji_qq * Pi;
    end

    % Arrange solution vector for ODE solvers
    qdot = zeros(max(size(qq),1));
    for i = 1:max(size(qq))
        qdot(i, 1) = qdoti(i);
    end
end