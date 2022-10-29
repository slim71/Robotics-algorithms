function points = BallNeighborhood(center, radius, n_samples)
    
    n_coord = length(center);
    points = center + radius .* rand(n_samples, n_coord) - radius/2;
end
