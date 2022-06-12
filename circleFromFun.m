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

end