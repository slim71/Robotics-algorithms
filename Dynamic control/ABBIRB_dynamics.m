%% Startup

% Launches the file that prepares robot and trajectory
if ~exist('abbirb', 'var')
    control_algorithms
end

%% Symbolic variables and parameters
syms mass [6 1]
syms x [6 1]
syms y [6 1]
syms z [6 1]
syms I [3 3 6]
syms q [6 1]
assume(q, 'real')
syms dq [6 1]
assume(dq, 'real')
syms ddq [6 1]
assume(ddq, 'real')
syms taus [6 1]

num_links = max(size(abbirb.links));
num_link_param = 10;

%% Dynamics

% Inertia
jac = abbirb.jacob0(q0);
Bdyn = zeros(num_links);
jj = zeros(num_links);
for i = 1:num_links
    jj(:, i) = jac(:, i);
    Bdyn = Bdyn + mass(i)*jj(1:3,:)'*jj(1:3,:) + jj(4:6,:)'*I(:,:,i)*jj(4:6,:);
end

% Coriolis
ch = @(i, j, k) 1/2 * (diff(Bdyn(i, j), q(k)) + diff(Bdyn(i, k), q(j)) - diff(Bdyn(j, k), q(i)));
Cdyn = zeros(num_links);
for i = 1:num_links
    for j = 1:num_links
        for k = 1:num_links
            Cdyn(i, j) = Cdyn(i,j) + ch(i,j,k)*dq(k);
        end
    end
end

% Gravity
gravity = [0, 0, -9.81]';
Gdyn = zeros(num_links, 1);
jj = zeros(num_links);
for i = 1:num_links
    jj(:, i) = jac(:, i);
    Gdyn = Gdyn - mass(i)*jj(1:3, :)' * gravity;
end

% Dynamics equation
dyn = Bdyn*ddq + Cdyn*dq + Gdyn;

%% Regressor

for i = 1:num_links
    parameters(i, 1) = mass(i);
    parameters(i, 2) = x(i);
    parameters(i, 3) = y(i);
    parameters(i, 4) = z(i);
    parameters(i, 5) = I(1,1, i);
    parameters(i, 6) = I(1,2, i);
    parameters(i, 7) = I(1,3, i);
    parameters(i, 8) = I(2,2, i);
    parameters(i, 9) = I(2,3, i);
    parameters(i, 10) = I(3,3, i);
end

for eq = 1:num_links
    for lin = 1:num_links
        for par = 1:num_link_param
            regY(eq, (lin-1)*num_link_param+par) = diff(dyn(lin), parameters(lin, par));
        end
    end
end