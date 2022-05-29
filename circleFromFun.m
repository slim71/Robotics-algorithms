function [xc, yc, zc] = circleFromFun(fun, time)

    Xs = [];
    Ys = [];
    Zs = [];

    for i = time
        Xs = [Xs; indexAt(fun(i), 1)];
        Ys = [Ys; indexAt(fun(i), 2)];
        Zs = [Zs; indexAt(fun(i), 3)];
    end

    xc = Xs;
    yc = Ys;
    zc = Zs;

%     figure()
%     hold on
%     pp = arrayfun(fun, time, 'UniformOutput', 'false');
%     p = cell2mat(pp);
%     plot3(p(:,1), p(:,2), p(:,3), 'ro-');
end