% take away cell cubes from total volume

function[ vol ] = compute_pareto_volume( cells, lower, upper )

vol = prod(upper-lower);

for s=1:length(cells)
    l  = cells{s}(1,:);
    u  = cells{s}(2,:);
    cell_vol = prod(u-l); 
    vol      = vol - cell_vol;
end



