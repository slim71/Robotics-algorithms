function dist = DlinearizedX0(f, G, vars, body_pars, init_cond, body_values)

	m = size(G, 2);
	n = size(f, 1);
	
	n_iter = 100;
	count = 0;
    % To keep track of previously computed brackets
	brackets = zeros(n, n_iter);

    % Starting off with input vectors gi already inside
    dist = subs(G, [vars, body_pars], [init_cond, body_values]);

    % Generating [f, gi], [f, [f, gi]], ....
	for i = 1:m
        % Must continue "indefinitly" if size keeps up
		while size(LinIndCols(brackets), 2) < n && (count < n_iter)
			count = count + 1;
			if count == 1 % First simple bracket
				brackets(:, count) = double(subs(LieBracket(f, G(:, count), vars), [vars, body_pars], [init_cond, body_values]));
            else % Composed brackets
				brackets(:, count) = double(subs(LieBracket(f, brackets(:, count-1), vars), [vars, body_pars], [init_cond, body_values]));
			end
		end
		dist = [dist brackets];
		brackets = zeros(n, n_iter);
		count = 0;
	end
	dist = LinIndCols(dist);
end
