function [ optimum ] = globalOptimization(target_noderiv, target, xmin, xmax)

	d = length(xmin);

	% We evaluate the objective in a grid and pick the best location

	gridSize = 100;

	Xgrid = repmat(xmin, gridSize, 1) + repmat(xmax - xmin, gridSize, 1) .* rand(gridSize, d);

    y = zeros(gridSize,1);
    for g=1:gridSize
    	y(g) = target_noderiv(Xgrid(g,:));
    end
    
	[ minValue minIndex ] = min(y);

	start = Xgrid(minIndex,:);

	% We optimize starting at the best location

	optimum = fmincon(target, start, [], [], [], [], xmin, xmax, [], ...
                optimset('MaxFunEvals', 50, 'TolX', 0.05, 'Display', 'off', 'GradObj', 'on'));

end
