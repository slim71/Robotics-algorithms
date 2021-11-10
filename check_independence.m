function ind = check_independence(a, b)
    together = [a, b];
    ind = ~(rank(together) < min(size(together)));
end