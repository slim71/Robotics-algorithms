function des_q = desiredConfiguration(qinit, pose, robot)
% Generate desired pose, given the initial configuration and

	n_joints = size(robot.links, 2);
    n_points = size(pose, 1);

    % Init parameters    
    des_q = zeros(n_joints, n_points);
    des_q(:, 1) = qinit;

    for i = 2:n_points
        pos = pose(i-1, 1:3);
        rotmat = EulToMat(pose(i-1, 4:6), "ZYX");
        hom_SE3 = SE3(rotmat, pos);
    
        des_q(:, i) = robot.ikine(hom_SE3, 'q0', des_q(:, i-1)');
    end

end
