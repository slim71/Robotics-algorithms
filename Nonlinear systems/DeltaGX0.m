function dg = DeltaGX0(G, vars, body_pars, init_cond, body_values)
    
    m = size(G, 2);

    dg = [];

    % Excluding :
    %  - (i,j)=(1,1) because [gi, gi] = 0
    %  - j<i because [gi, gj] = -[gj, gi] --> not linearly indipendent
    for i = 1:m
        for j = i+1:m
            bracket = LieBracket(...G(:, i), G(:,j), vars);
                        double(subs(G(:, i), [vars, body_pars], [init_cond, body_values])), ...
                        double(subs(G(:, j), [vars, body_pars], [init_cond, body_values])), ...
                        vars);
            dg = LinIndCols([dg bracket]);
        end
    end
end
