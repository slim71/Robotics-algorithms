clc; clear all; close all;

syms t
syms q [6 1]

A     = 35/1000;
B     = 25/1000;
CC    = 80/1000;
F     = 515/1000;
H     = 400/1000;
G     = 560/1000;
Omega = 1000/1000;
Rho   = 150/1000;
d1    = -Rho;
d2    = B + F + CC;
d3    = H + G + A;

% q = zeros(6); %q(t)

denavit = [[B,  pi/2,   H,  pi/2 + q1,   "R"] % all qx(t)
           [G,  0,      0,  pi/2 + q2,   "R"]
           [A,  pi/2,   0,         q3,   "R"]
           [0,  -pi/2,  F,         q4,   "R"]
           [0,  pi/2,   0,         q5,   "R"]
           [0,  0,      CC,        q6,   "R"]];
       
Jac = simplify(JacobFromDH(denavit)); % Jac(t)
% Jac = symfun(simplify(JacobFromDH(denavit)), [q1 q2 q3 q4 q5 q6]);

% Just to have an example
% Commented out because too computationaly heavy
% pj = pinv(Jac);
% dj = damped(Jac, 100);

twistEE(t) = Omega*Rho*[cos(Omega*t); -sin(Omega*t); 0; 0; 0; 0]; 
Oee(t) = [d1+Rho*cos(Omega*t), d2, d3+sin(Omega*t)];
OeeHat(t) = AxisToSkew(Oee(t));
Moee(t) = [[eye(3),     OeeHat(t)];
           [zeros(3,3), eye(3)   ]];
twistO(t) = Moee(t)*twistEE(t);

qSotNorm = sot(Jac, twistO(t), 0);
% qSotNorm = sot(Jac(t), xDotDes(t), 0);
% qSotDamp = sot(Jac, xDotDes, 1, 100);

homDes(t) = HomX(0,[d1,d2,d3])*...
            HomX(0,[Rho*cos(Omega*t), 0, Rho*sin(Omega*t)])*...
            HomY(pi/2,[0,0,0])*...
            HomZ(pi/2,[0,0,0]);

qq = qSotNorm;

%% Time integration
tfinal = 2*pi/Omega;

%% Robot design
kuka = rigidBodyTree('Dataformat', 'column');%, 'MaxNumBodies', 6);

dhparams = [B  pi/2     H   pi/2+qq(1,2,1);
            G  0        0   pi/2+qq(2,2,1);
            A  pi/2     0   qq(3,2,1);
            0  -pi/2    F   qq(4,2,1);
            0  pi/2     0   qq(5,2,1);
            0  0        CC  qq(6,2,1)];

bodyNames = {'b1','b2','b3','b4','b5', 'b6'};
parentNames = {'base','b1','b2','b3','b4', 'b5'};
jointNames = {'j1','j2','j3','j4','j5', 'j6'};
jointTypes = {'revolute','revolute','revolute','revolute','revolute', 'revolute'};

for k = 1:6
    % Create a rigidBody object with a unique name
    kukaBodies(k) = rigidBody(bodyNames{k});
    % Create a rigidBodyJoint object and give it a unique name
    kukaBodies(k).Joint = rigidBodyJoint(jointNames{k}, jointTypes{k});
    % Use setFixedTransform to specify the body-to-body transformation using DH parameters
    setFixedTransform(kukaBodies(k).Joint, dhparams(k,:), 'dh');
    % Attach the body joint to the robot
    addBody(kuka, kukaBodies(k), parentNames{k});
end

% Add a final body to function as the end-effector (handle)
% bn = 'handle';
% ee = rigidBody(bn);
% setFixedTransform(ee.Joint,trvec2tform([0 -0.15 0]));
% addBody(kuka, ee, 'b5');

showdetails(kuka)

f = figure2();
ax = show(kuka, qq(:, end, 1));
p = uipanel(f, 'Units', 'normalized','Position',[0.7 0.04 0.25 0.05]);
c = uicontrol(p,'Style','slider', 'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.7]);
c.Value = 1;
c.Min = 1;
c.Max = size(qq,3);
c.SliderStep = [1/c.Max 10/c.Max];
addlistener(c, 'Value', 'PostSet', @(hObject, events) animation_sot(hObject, events, kuka, qq));


qq1 = mod(squeeze(qq(1,2,:)),2*pi);
qq2 = mod(squeeze(qq(2,2,:)),2*pi);
qq3 = mod(squeeze(qq(3,2,:)),2*pi);
qq4 = mod(squeeze(qq(4,2,:)),2*pi);
qq5 = mod(squeeze(qq(5,2,:)),2*pi);
qq6 = mod(squeeze(qq(6,2,:)),2*pi);
f2 = figure2();
hold on
grid
plot(qq1)
plot(qq2)
plot(qq3)
plot(qq4)
plot(qq5)
plot(qq6)

% 

% hold on
% tt = 0:1:10;
% plot3(zeros(1,11),Omega*Rho*cos(Omega*tt), Omega*Rho*(-sin(Omega*tt)))
% circle = [zeros(1,11);Omega*Rho*cos(Omega*tt); Omega*Rho*(-sin(Omega*tt));ones(1,11)]
% app = mat*circle
% plot3(app(1,:), app(2,:), app(3,:))
% hold on
% plot3(app(1,:), app(2,:), app(3,:))