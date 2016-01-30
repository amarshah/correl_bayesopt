% Y is   NxR

function[ cells, pareto ] = obs_to_integration_cells( Y, lower, upper )

% compute pareto front
pareto = pareto_frontier( Y );     % Np x R

if length(pareto)== 0
    pareto
end

% include upper and lower limits
points = [lower;pareto;upper];

% compute all the grid cells C 
C = all_grid_cells( points );

% keep relevant grid cells
cells = keep_relevant_cells( C, pareto );






