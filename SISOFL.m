function [r, Lfs, Lgs] = SISOFL(f, g, h, vars, init_vars)

    % Init
    r = 0;

    % This are meant to be taken as 
    %   (LF^i)h     , i=0...r
    %   (LgLf^(i))h , i=0...r-1
    Lfs(1) = h;
    Lgs = sym.empty;
    
    while r < length(f)
        r = r + 1;

        lf_i = ScalarFDerivative(f, Lfs(r), vars);
        lglf_i = ScalarFDerivative(g, Lfs(r), vars);

        % Append iteration results and continue
        Lfs(end+1) = lf_i;
        Lgs(end+1) = lglf_i;

        % If LgLf h(x_0) != 0, r is the relative degree
        if lglf_i ~= 0 %subs(lglf_i, vars, init_vars) ~= 0
            break
        end

    end

    if subs(lglf_i, vars, init_vars) == 0
        % Relative degree not defined
        r = r * (-1);
    end

end
