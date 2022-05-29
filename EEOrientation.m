function EEOrientation(point, directions)
% TODO

    % Get axis directions
    x_dir = directions(:, 1);
    y_dir = directions(:, 2);
    z_dir = directions(:, 3);

    % Plot quivers along axis
    quiv_x = quiver3(point(1), point(2), point(3), x_dir(1), x_dir(2), x_dir(3));
    quiv_x.Color = 'red';
    quiv_y = quiver3(point(1), point(2), point(3), y_dir(1), y_dir(2), y_dir(3));
    quiv_y.Color = 'green';
    quiv_z = quiver3(point(1), point(2), point(3), z_dir(1), z_dir(2), z_dir(3));
    quiv_z.Color = 'blue';

end