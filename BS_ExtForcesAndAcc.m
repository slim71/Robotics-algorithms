function [torque, accel] = BS_ExtForcesAndAcc(robot, q, dq, desq, desdq, desddq, lambda, kp, kd)
    
    % Ensure input arrays have proper dimensions
    assert(isrow(q) && isrow(dq) && isrow(desq) && isrow(desdq) && isrow(desddq), ...
        "All input arrays must be rows!");

    % Errors
    e = desq - q;
    de = desdq - dq;

    dq_ref = desdq' + lambda*e';
    ddq_ref = desddq' + lambda*de';
    s = de' + lambda*e';

    M = robot.inertia(q);
    C = robot.coriolis(q, dq);
    G = robot.gravload(q)';

    torque =  M*ddq_ref + C * dq_ref + G + lambda*s + e';
 
    accel = pinv(M) * (torque - C*dq' - G);

end
