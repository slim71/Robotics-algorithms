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

xDotDes(t) = Omega*Rho*[cos(Omega*t); -sin(Omega*t); 0; 0; 0; 0]; 

qSotNorm = sot(Jac, xDotDes(t), 0);
% qSotNorm = sot(Jac(t), xDotDes(t), 0);
% qSotDamp = sot(Jac, xDotDes, 1, 100);

qq = qSotNorm;

%% Time integration
tfinal = 2*pi/Omega;

%% Robot design
kuka = rigidBodyTree('Dataformat', 'column', 'MaxNumBodies', 6);

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

    kukabody = rigidBody(bodyNames{k});
    kukabody.Joint = rigidBodyJoint(jointNames{k}, jointTypes{k});
    
    setFixedTransform(kukabody.Joint, dhparams(k,:), 'dh');
    
    if k == 1
        mat = [[1,  0,  0,  -3/20 + 3*cos(2*pi)/20]
               [0,  0, -1,  31/50]
               [0,  1,  0,  199/200 + 3*sin(2*pi)/20]
               [0,  0,  0,  1]];
        
       setFixedTransform(kukabody.Joint, mat*kukabody.Joint.JointToParentTransform)
    end
    
    
    
    addBody(kuka, kukabody, parentNames{k});
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

% 

% hold on
% tt = 0:1:10;
% plot3(zeros(1,11),Omega*Rho*cos(Omega*tt), Omega*Rho*(-sin(Omega*tt)))
% circle = [zeros(1,11);Omega*Rho*cos(Omega*tt); Omega*Rho*(-sin(Omega*tt));ones(1,11)]
% app = mat*circle
% plot3(app(1,:), app(2,:), app(3,:))
% hold on
% plot3(app(1,:), app(2,:), app(3,:))