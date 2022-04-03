%% Clik kinematics

[t_clik, y_clik] = ode15s(@(t,y) clik(t,y,J,destwist), [t0 tf], qinit);
q_clik = y_clik';

ode_fig = figure2('Name', 'Clik with ode45');
for i=[1:size(q_clik, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_clik(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
clik_fig = figure2('Name', 'Clik solution');
for j=[1:size(q_clik,1)]
    subplot(size(q_clik,1), 1, j)
    plot(q_clik(j,:))
end

%% Clik errors
pe_fig = figure2('Name', 'Clik position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_clik,1)]
    pe = pos_error(t_clik(i), q_clik(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_clik,1)]   
    subplot(3, 1, 1);
    plot(t_clik, p1);
    subplot(3, 1, 2);
    plot(t_clik, p2);
    subplot(3, 1, 3);
    plot(t_clik, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'Clik orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_clik,1)]
    oe = or_error(t_clik(i), q_clik(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_clik,1)]   
    subplot(3, 1, 1);
    plot(t_clik, o1);
    subplot(3, 1, 2);
    plot(t_clik, o2);
    subplot(3, 1, 3);
    plot(t_clik, o3);
end

%% Clik with errors
gains_cwe = 50 * [1 1 1 1 1 1];

[t_cwe, y_cwe] = ode15s(@(t,y) clik_with_error(t,y,J,destwist,tot_errors, gains_cwe),...
                       [t0 tf], qinit);
q_cwe = y_cwe';

ode_fig2 = figure2('Name', 'Cwe with ode45');
for i=[1:size(q_cwe, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_cwe(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
cwe_fig = figure2('Name', 'CWE solution');
for j=[1:size(q_cwe,1)]
    subplot(size(q_cwe,1), 1, j)
    plot(q_cwe(j,:))
end

%% Clik_with_errors errors
pe_fig = figure2('Name', 'CWE position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_cwe,1)]
    pe = pos_error(t_cwe(i), q_cwe(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_cwe,1)]   
    subplot(3, 1, 1);
    plot(t_cwe, p1);
    subplot(3, 1, 2);
    plot(t_cwe, p2);
    subplot(3, 1, 3);
    plot(t_cwe, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'CWE orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_cwe,1)]
    oe = or_error(t_cwe(i), q_cwe(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_cwe,1)]   
    subplot(3, 1, 1);
    plot(t_cwe, o1);
    subplot(3, 1, 2);
    plot(t_cwe, o2);
    subplot(3, 1, 3);
    plot(t_cwe, o3);
end
%% Clik with secondary objective
gains_cso = 50 * [1 1 1 1 1 1];

k0 = eye(7);
mins = [-pi, -pi, -pi, -pi, -pi, -pi, -pi];
maxs = [pi, pi, pi, pi, pi, pi, pi];

H = @(q) JointCenterRange(q, mins, maxs);

[t_cso, y_cso] = ode45(@(t,y) clik_secobj(t,y,J,destwist,tot_errors, gains_cso, H),...
                       [t0 tf], qinit);
q_cso = y_cso';

ode_fig2 = figure2('Name', 'Cso with ode45');
for i=[1:size(q_cso, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_cso(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
cwe_fig = figure2('Name', 'CSO solution');
for j=[1:size(q_cso,1)]
    subplot(size(q_cso,1), 1, j)
    plot(q_cso(j,:))
end

%% Clik_secobj errors
pe_fig = figure2('Name', 'CSO position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_cso,1)]
    pe = pos_error(t_cso(i), q_cso(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_cso,1)]   
    subplot(3, 1, 1);
    plot(t_cso, p1);
    subplot(3, 1, 2);
    plot(t_cso, p2);
    subplot(3, 1, 3);
    plot(t_cso, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'CSO Orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_cso,1)]
    oe = or_error(t_cso(i), q_cso(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_cso,1)]   
    subplot(3, 1, 1);
    plot(t_cso, o1);
    subplot(3, 1, 2);
    plot(t_cso, o2);
    subplot(3, 1, 3);
    plot(t_cso, o3);
end

%% Clik with q halfway
gains_half = 50 * [1 1 1 1 1 1];
mins = [-pi, -pi, -pi, -pi, -pi, -pi, -pi];
maxs = [pi, pi, pi, pi, pi, pi, pi];

[t_half, y_half] = ode45(@(t,y) clik_secobj(t,y,J,destwist,tot_errors, gains_half, mins, maxs),...
                       [t0 tf], qinit);
q_half = y_half';

ode_half = figure2('Name', 'Halfq with ode');
for i=[1:size(q_half, 2)]
    plot3(first_traj(:,1),first_traj(:,2),first_traj(:,3),'-', 'Color', 'black')
    hold on
    grid
    show(kuka, q_half(:,i));
    view(45,45)
    drawnow
    pause(0.1)
    hold off
end
half_fig = figure2('Name', 'CSO solution');
for j=[1:size(q_half,1)]
    subplot(size(q_half,1), 1, j)
    plot(q_half(j,:))
end

%% Clik_halfq errors
pe_fig = figure2('Name', 'CSO position error');
p1 = [];
p2 = [];
p3 = [];
for i=[1:size(t_half,1)]
    pe = pos_error(t_half(i), q_half(:,i));
    p1 = [p1, pe(1)];
    p2 = [p2, pe(2)];
    p3 = [p3, pe(3)];
end

for i=[1:size(t_half,1)]   
    subplot(3, 1, 1);
    plot(t_half, p1);
    subplot(3, 1, 2);
    plot(t_half, p2);
    subplot(3, 1, 3);
    plot(t_half, p3);
end

% Orientation error computed using quaternions

or_fig = figure2('Name', 'Halfq Orientation error');
o1 = [];
o2 = [];
o3 = [];
for i=[1:size(t_half,1)]
    oe = or_error(t_half(i), q_half(:,i));
    o1 = [o1, oe(1)];
    o2 = [o2, oe(2)];
    o3 = [o3, oe(3)];
end

for i=[1:size(t_half,1)]   
    subplot(3, 1, 1);
    plot(t_half, o1);
    subplot(3, 1, 2);
    plot(t_half, o2);
    subplot(3, 1, 3);
    plot(t_half, o3);
end

