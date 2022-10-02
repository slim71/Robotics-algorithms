function [rs, Lfs, Gamma, E] = MIMOFL(f, gs, hs, vars, init_vars)

    % Init
    rs = zeros(size(hs,1), size(gs, 2));

    % This are meant to be taken as 
    %   (LF^i)h     , i=0...r
    %   (LgLf^(i))h , i=0...r-1
    Gamma = sym.empty;
    E = sym.empty;
    Lfs = cell(1, size(hs, 1));
    
    for i = 1:size(hs, 1)

        % ris: 1xm
        % Lfi: 1xr_i
        % Lgi: mxr_i
        [ris, Lfi, Lgi] = MISOFL(f, gs, hs(i), vars, init_vars);

        rs(i, :) = ris;

        Lfs(i) = {Lfi(1:end-1)}; % Return empty cell if Lfi is empty
        if isempty(Lfi)
            Gamma(i, 1) = 0;
        else
            Gamma(i, 1) = Lfi(end);
        end
        E(i, :) = Lgi(:, end)';

    end
end