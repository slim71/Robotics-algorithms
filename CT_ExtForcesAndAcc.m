function [torque, accel] = CT_ExtForcesAndAcc(robot, q, dq, desq, desdq, desddq, Kp, Kd)
    
    % Ensure input arrays have proper dimensions
    assert(isrow(q) && isrow(dq) && isrow(desq) && isrow(desdq) && isrow(desddq), ...
        "All input arrays must be rows!");

    % Errors
    err = (desq - q);
    derr = (desdq - dq);

    % Robot dynamics matrices
    M = robot.inertia(q);
    C = robot.coriolis(q, dq);
    G = robot.gravload(q)';

    torque =  M * (desddq' + Kd*derr' + Kp*err') + C * dq' + G;
    accel = pinv(M) * (torque - C*dq' - G);

end
