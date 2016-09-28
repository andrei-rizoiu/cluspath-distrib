function legends = plot_graphic_2measures(History, xmeasure, ymeasure, bestPoint, singleGraphic)

	if (nargin < 5)
		singleGraphic = true;
	endif
	
	baseMarkerSize = 2;
	baseLineWidth = 2;
	if (singleGraphic)
		baseMarkerSize = 4;
		baseLineWidth = 4;
	endif

	colors = {"blue", "black", "cyan", "green", "magenta", "red", "yellow"};

	yval = [];
	xval = [];
	paretopoints = [];
	legends = {};
	survived = [];
	colindex = 0;
	for i=1:(length(History)-1)
		colindex += 1;
		if (colindex > 7)
			colindex -= 7;
		endif
		
		xval = History(i).Measures(xmeasure, :) ;
		yval = History(i).Measures(ymeasure, :) ;
		survived = cell2mat(struct2cell(History(i).Outputs)(5,1,:)(:));
		paretopoints = [ paretopoints ; History(i).Measures([xmeasure ymeasure], survived == 1)'] ;
		legends = [ legends {sprintf("Generation %d", i)}];
		
		plot(xval(survived == 1), yval(survived == 1), '+', 'LineWidth', baseLineWidth, 'markersize', baseMarkerSize, 'color', colors{colindex});
		hold on;
		plot(xval(survived == 0), yval(survived == 0), 'o', 'LineWidth', baseLineWidth, 'markersize', baseMarkerSize, 'color', colors{colindex});
	endfor
	
	plot(History(end).Measures(xmeasure, :), History(end).Measures(ymeasure, :), '^r', 'LineWidth', baseLineWidth+1, 'markersize', baseMarkerSize+1);
	paretopoints = [ paretopoints ; History(end).Measures([xmeasure ymeasure], :)'] ;
	legends = [ legends {"Final generation"}];
	
	# here, I can use the getParetoPoints function, which is designed for 2-d pareto fronts.
	[pareto, paretoindex] = getParetoFront(paretopoints);
	plot(pareto(:, 1), pareto(:, 2), "-b");
	legends = [ legends {"Pareto front"}];
	
	if (nargin >= 4 && ~isempty(bestPoint))
		# best point was given to be ploted
		plot(bestPoint(xmeasure), bestPoint(ymeasure), '*b', 'LineWidth', baseLineWidth+2, 'markersize', baseMarkerSize+3); 
	endif
	
	indics = {"MDvar" "Tvar" "ShaP" "SPass"};
	xlabel(indics{xmeasure});
	ylabel(indics{ymeasure});
	title(sprintf("Evolutionary populations, %s-%s measures", indics{xmeasure}, indics{ymeasure}));
	
	grid on;
endfunction
