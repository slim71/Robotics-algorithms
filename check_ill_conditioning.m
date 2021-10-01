function [ill, num_cond] = check_ill_conditioning(mat)
% threshold at 10^3

    num_cond = cond(mat);
    ill = (num_cond >= 10^3);
end