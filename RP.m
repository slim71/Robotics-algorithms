function qsol = RP(J, xdot, kind)
    qdot(l+1) = 0;
    P(l+1) = eye();
    
    % J not computable, so neither P is
    if (kind == 0)
        for task = [l:-1:1]
            P(task+1) = I - pinv(J(task)) * J(task);
            qdot(task) = qdot(task-1) + pinv(J(task) * P(task+1)) * (xdot(task)-J(task) * qdot(task+1));
        end
    else
        for task = [l:-1:1]
            P(task+1) = I - damped(J(task)) * J(task);
            qdot(task) = qdot(task-1) + damped(J(task) * P(task+1)) * (xdot(task)-J(task) * qdot(task+1));
        end
    end
    
    qsol = qdot(end);
end