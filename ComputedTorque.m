function [torque, accel] = ComputedTorque(robot, q, dq, desq, desdq, desddq)
    
    % Ensure input arrays have proper dimensions
    assert(isrow(q) && isrow(dq) && isrow(desq) && isrow(desdq) && isrow(desddq), ...
        "All input arrays must be rows!");

    % Compute torque and acceleration
    Kp = 20 * diag([1 1 1 1 1 1]);
    Kd = 10 * diag([1 1 1 1 1 1]);

    err = (desq - q);
    derr = (desdq - dq);

    M = robot.inertia(q);
    C = robot.coriolis(q, dq);
    G = robot.gravload(q)';

    torque =  M * (desddq' + Kd*derr' + Kp*err') + C * dq' + G;
 
    accel = pinv(M) * (torque - C*dq' - G);

end
