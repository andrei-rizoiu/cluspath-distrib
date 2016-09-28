function [bestPoint, bestParameters, paretoplot] = determine_best_individual(varargin)

	[History, optimizationCriteria, generation ] = ... 
		process_options (varargin, 'History', [], 'optimizationCriteria', [1 2 3 4], 'generation', [] );

	if ( isempty(History))
		error("I cannot work with an empty History. Give me the data structure or the file in which it is saved.");
	endif
	
	# if here, we were given a History to continue
	# is it a structure (an Octave object?)
	if (isstruct(History))
		ng = length(History);
		Configs = History(ng).Configs;
		Outputs = History(ng).Outputs;
		Measures = History(ng).Measures;
	else
		# maybe it is a file to load from; lets see what is in it?
		[ vars ] = who("-file", History);
		correct = length(vars) == 1;	%need to have the 1 variables here after
		for i=1:length(vars)
			if ( strcmp( vars{i}, "History" ) == 0 )
				correct = false;
			endif
		endfor
	
		if ( correct )
			load( History );
			ng = length(History);
			Configs = History(ng).Configs;
			Outputs = History(ng).Outputs;
			Measures = History(ng).Measures;
		else
			error(sprintf("The given file '%s' is not a valid history file!"));
		endif
	endif
	
	if (isempty(generation))
#		generation = [1:length(History)];
		generation = length(History);
	endif
	
	############################ WORK #################################
	# find out best point (point on the Pareto front, closest to the origin)
	bestDist = inf;
	bestGeneration = 0;
	bestPoint = [];
	points = [];
	for i=generation
		points = [ points ...
				[	History(i).Measures(optimizationCriteria, :); ...
					repmat([i], 1, columns(History(i).Measures(optimizationCriteria, :))); ...
					[1:columns(History(i).Measures(optimizationCriteria, :))] ... 
				] ...
			] ;
	endfor	

	## determine the Pareto front - cannot use the getParetoFront function, since it works only in 2D
#	[paretoplot, paretoindex] = getParetoFront(points(1:end-2, :)');
	counts = domination_count(points(1:end-2, :), [1:length(optimizationCriteria)]);
	points = points(:, counts == 0);
	paretoplot = points(1:end-2, counts==0)';
	vals = normalizeMatr(points(1:end-2, :)', false)'; #normalize the values so that the scale differences don't influence too much. No standardize
	
	center = zeros(rows(vals), 1);
	my_dist = distance(vals, center);
	[bestDist, pos] = min(my_dist);
	bestGeneration = points(end-1, pos);
	bestPoint = History(bestGeneration).Measures(:, points(end, pos));
	bestConfig = History(bestGeneration).Configs(points(end, pos)).params;
	
	printf("**> Best overall individual was in generation %d/%d: \n", bestGeneration, length(History));
	printf("\t\t\tOptimizing following criteria: [");
	if (ismember(1, optimizationCriteria))
		printf("MDvar");
	endif
	if (ismember(2, optimizationCriteria))
		printf(", Tvar");
	endif
	if (ismember(3, optimizationCriteria))
		printf(", ShaP");
	endif
	if (ismember(4, optimizationCriteria))
		printf(", SPass");
	endif
	printf("]\n");
	printf("\t\t\tMeasures:\tMDvar=%.3f, Tvar=%.3f, ShaP=%.3f, SPass=%.3f\n", bestPoint(1), bestPoint(2), bestPoint(3), bestPoint(4)); 
	printf("\t\t\tParameters:\tAlpha=%.3f, Beta=%.4f, Delta=%.3f, lambda1=%.2f, lambda2=%.2f, lambda3=%.2f\n", ...
	    bestConfig.Alpha, bestConfig.Beta, bestConfig.Delta, bestConfig.lamb1, bestConfig.lamb2, bestConfig.lamb3);
    		
	bestParameters = [bestConfig.Alpha, bestConfig.Beta, bestConfig.Delta, bestConfig.lamb1, bestConfig.lamb2, bestConfig.lamb3];
    	
endfunction
