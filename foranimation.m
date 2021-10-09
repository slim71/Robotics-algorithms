[X,Y] = meshgrid(linspace(-5, 5, 50));
fcn = @(x,y,k) k*x.^2 + y.^2;
v = [1:-0.05:-1;  -1:0.05:1];
for k1 = 1:2
    for k2 = v(k1,:)
        surfc(X, Y, fcn(X,Y,k2))
        axis([-5  5    -5  5    -30  50])
        drawnow
        pause(0.1)
    end
end

%% Create Animation of Streaming Data
% Create an animation of a line growing as it accumulates 2,000 data points. 
% Use |drawnow| to display the changes on the screen after each iteration through 
% the loop.

h = animatedline;
axis([0 4*pi -1 1])
x = linspace(0,4*pi,2000);

for k = 1:length(x)
    y = sin(x(k));
    addpoints(h,x(k),y);
    drawnow
end
%% 
% Copyright 2015-2017 The MathWorks, Inc.
