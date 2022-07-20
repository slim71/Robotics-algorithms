function [torque, accel] = BackStepping(robot, q, dq, desq, desdq, desddq, gain)
    
    % Ensure input arrays have proper dimensions
    assert(isrow(q) && isrow(dq) && isrow(desq) && isrow(desdq) && isrow(desddq), ...
        "All input arrays must be rows!");

    % Errors
    err = (desq - q);
    s = (desdq - dq);

    M = robot.inertia(q);
    C = robot.coriolis(q, dq);
    G = robot.gravload(q)';
    J = robot.jacob0(q);

    torque =  M*desddq' + C * desdq' + G + gain*s' + J'*err';
 
    accel = pinv(M) * (torque - C*dq' - G);

end
