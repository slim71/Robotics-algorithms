function [flag, cond] = STLCCheck(f, G, vars, body_pars, init_cond, body_values)
% STLCCHECK
%
% The following are sufficient conditions to affirm that a system which
%	is s.t.l.a. at x_0 is also s.t.l.c. (n: system dim, m: input dim):
%		i)	 f(x)=0 for all x inside B(x_0)
%		ii)	 f(x) inside span{gi(x),..., gm(x)} for all x inside B(x_0)
%		iii) f(x_0)=0 
%               AND
%            dim(Dlin(x0))=dim(span{gi,[f, gi], [f,[f, gi]], ..., [(ad^mu_i)_f](g), ...})=n
%		        for each i = 1...m and mu_i in N 
%               -> this is because Dlinearized ~ [B, AB, ..., A^nB] + O(|x^2|)
%		iv)  f(x_0)=0
%               AND
%            dim(DL(x_0)+DG(x_0))=n with DG=span{[gi, gj]} 
%				for each i,j = 1...m
%		v)	 f(x_0)=0
%               AND
%            every "bad brackets" (brackets with odd # f and even # g, not null)
%            is a linear combination of "good brackets" (smaller total # f+g)
%		vi)	...

    n = size(f, 1);
    m = size(G, 2);

    stlc = false;
    condition = 0;

    % Neighborhood of x0, made from randomly selected points
    hood_dim = 100;
    hood_radius = 0.01;
    thehood = BallNeighborhood(init_cond, hood_radius, hood_dim);

    % Base: f(x_0) = 0
    f_x0 = subs(f, [vars, body_pars], [init_cond, body_values]);
    if f_x0 ~= 0
        flag = false;
        return
    end

    while ~stlc && condition <= 6
        switch(condition)

            case 1
                % First condition: 
                %   f(x)=0 for all x inside B(x_0)
                stlc = true;
                for i = 1:size(hood_dim)
                    f_xi = subs(f, [vars, body_pars], [thehood(i, :), body_values]);
                    if double(f_xi) ~= 0
                        stlc = false;
                        break
                    end
                end

            case 2
                % Second condition: 
                %   f(x) inside span{gi(x),..., gm(x)} for all x inside B(x_0)
                stlc = true;
                for i = 1:size(hood_dim)
                    f_xi = subs(f, [vars, body_pars], [thehood(i, :), body_values]);
                    if rank([G f_xi]) > rank(G)
                        stlc = false;
                        break
                    end
                end

            case 3
                % Third condition
                %   f(x_0)=0 
                %      AND
                %   dim(Dlin(x0)) = dim(span{gi,[f, gi], [f,[f, gi]], ..., [(ad^mu_i)_f](g), ...})=n
	            %        for each i = 1...m and mu_i in N (n: system dim, m: input dim)
                %      -> this is because Dlinearized ~ [B, AB, ..., A^nB] + O(|x^2|)
                stlc = true;
                DL = DlinearizedX0(f, G, vars, body_pars, init_cond, body_values);
                dl_dim = size(DL, 1);
                if dl_dim ~= n
                    stlc = false;
                end

            case 4
                % Fourth condition
                %   f(x_0)=0 and dim(DL(x_0)+DG(x_0))=n with DG=span{[gi, gj]} 
				%       for each i,j = 1...m
                stlc = true;
                DG = DeltaGX0(G, vars, body_pars, init_cond, body_values);
                dg_dim = size (DG, 2);
                if dl_dim + dg_dim ~= n
                    stlc = false;
                end

            case 5 % TODO: not fully tested! BEWARE
                % Fifth condition
                % f(x_0)=0
                %   AND
                % every "bad brackets" (brackets with odd # f and even # g, not null)
                % is a linear combination of "good brackets" (smaller total # f+g)
                stlc = true;
                
                deep_iters = 10;

                % To keep track of previously computed brackets
                % Starting from [f, G]
	            norm_brackets(:, :, 1) = double(subs(LieBracketsDistribution(f, G, vars), ...
                                [vars, body_pars], [init_cond, body_values]));
                bad_brackets = double.empty(size(norm_brackets, 1), 100, 0);

                i = 0;
            
                % Generating [f, G], [f, [f, G]], [G, [f, G]], ....
                while i < deep_iters
                    i = i + 1;

                    if mod(i, 2) == 0 % i is even --> no bad brackets will be generated

                        % Select previous step brackets
                        % i/2 normal and i/2 bad brackets were generated at last loop step
                        prev_norm_brackets = norm_brackets(:, :, end-(i/2-1):end);
                        prev_bad_brackets = bad_brackets(:, :, end-(i/2-1):end);

                        % Generate new brackets
                        for j = 1:i/2
                            % Every one from the previous normal brackets
                            norm_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(f, prev_norm_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                            norm_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(G, prev_norm_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                
                            % this is unique and not generated only for the
                            % last bad bracket, the one with higher # of g
                            if j == i/2
                                norm_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(G, prev_bad_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                            end
                        end
                        
                    else % i is odd
                        % all normal brackets were generated at last loop step
                        prev_loop_brackets = norm_brackets(:, :, end-(i-1):end);
                        current_bad_brackets = double.empty(size(prev_loop_brackets, 1),100,0);

                        % Create subspace to check into before adding new
                        % brackets
                        norm_subspace = reshape(norm_brackets, size(norm_brackets, 1), ...
                                           size(norm_brackets, 3)*size(norm_brackets,2));
                        bad_subspace = reshape(bad_brackets, size(bad_brackets, 1), ...
                                           size(bad_brackets, 3)*size(bad_brackets,2));
                        subspace = [norm_subspace bad_subspace];

                        % Generate new brackets
                        for j = 1:i
                            % this is unique and to be generated only for
                            % the first previous bracket, the one with
                            % higher # of f
                            if j == 1 
                                norm_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(f, prev_loop_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                            end

                            if mod(j, 2) == 0 % j is even
                                norm_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(G, prev_loop_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                            else % j is odd
                                % Do not add right away: check if condition
                                % is met first
                                current_bad_brackets(:, :, end+1) = double(subs(LieBracketsDistribution(G, prev_loop_brackets(:,:,j), vars), [vars, body_pars], [init_cond, body_values]));
                            end
                        end

                        % Check if bad brackets are linearly dependent from
                        % smaller good brackets
                        prev_rank = rank(subspace);

                        % Using the number of bad brackets generated right now
                        for j = 1:size(current_bad_brackets, 3)
                            % True if bad brackets is non linearly
                            % dependent from good brackets
                            if rank([subspace current_bad_brackets(:, :, j)]) > prev_rank
                                stlc = false;
                                break;
                            else
                                % add newly created bad brackets to the list
                                bad_brackets(:, :, end+1) = current_bad_brackets(:, :, j);
                            end
                        end

                    end

                end

            case 6 % TODO
        end

        condition = condition + 1;
    end

    flag = stlc;
    cond = condition - 1;

end
