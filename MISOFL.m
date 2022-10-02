function [rin_mat, LF_in, LGs_in] = MISOFL(f, gs, h, vars, init_vars)

    % Init
    rin_mat = zeros(1, size(gs, 2));
    min_rin = length(f);
    min_abs_r = min_rin;

    % This are meant to be taken as 
    %   (LF^i)h     , i=0...r
    %   (LgLf^(i))h , i=0...r-1
    LF_in = sym.empty;
    LGs_in = sym.empty;
    
    % Call SISO function for each input
    for j = 1:size(gs, 2)

        [rin_mat(j), Lfj, Lgj] = SISOFL(f, gs(:, j), h, vars, init_vars);
        rj = rin_mat(j);

        % Definite relative degree
        if (rj >= 0) && (rj < min_rin)
            min_rin = rj;
            LF_in = Lfj;

            min_abs_r = min(min_abs_r, abs(rj));

            if j > 1
                LGs_in = LGs_in(:, 1:rj);
            end

        % Non-definite relative degree
        elseif (rj < 0) && (abs(rj) < min_abs_r)
            min_abs_r = abs(rj);

            if j > 1
                LGs_in = LGs_in(:, 1:abs(rj));
            end

        end

        if j == 1
            LGs_in = Lgj;
        else
            LGs_in(j, :) = Lgj(1:min_abs_r);
        end

    end

end
