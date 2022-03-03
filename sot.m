function qdot = sot(tt, qq, sym_Js, sym_twistes, lambda)
    % initialization
    q0dot = zeros(max(size(qq)),1);

    % init
    qdoti = q0dot;
    Pi = eye(max(size(qq)));
    tot_J = [];

    % i=2...n:
    % q'_i = q'_(i-1) + inv(J_i*P_(i-1)) * (task_i - J_i*q'_(i-1))
    % P_i = P_(i-1) - inv(J_i * P_(i-1)) * J_i * P_(i-1)
    for i=[1:max(size(sym_twistes))]
        Ji = sym_Js{i}; % function handle
        twisti = sym_twistes{i}; % function handle

        Ji_qq = Ji(qq);
        tot_J = [tot_J; Ji_qq];

        [log, cnum] = check_ill_conditioning(Ji_qq * Pi);
        if log
            inverted = damped_pinv(Ji_qq*Pi, lambda);
        else
            inverted = pinv(Ji_qq*Pi);
        end

        qdoti = qdoti + inverted * (twisti(tt)' - Ji_qq * qdoti);
        % P not so much different
        Pi = eye(max(size(qq))) - pinv(tot_J)*tot_J; %Pi - pinv(Ji_qq*Pi) * Ji_qq * Pi;
    end

    for i=[1:max(size(qq))]
        qdot(i,1) = qdoti(i);
    end
end