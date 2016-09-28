## Usage:
##
## Function which attempts to find the best values for the parameters of the temporal clustering, using evolutionary algorithms.
## Each individual in the population is a solution, denoted by a set of parameters and the centroids themselves. The determining
## of the parameters is done by constructing the 4-dimentional Pareto front. Each individual contains 6 parameters of the algoritm:
## (Alpha, Beta, delta, lambda1, lambda2, lambda3) and the set of centroids issued after a couple (number of iterations - configurable)
## of clustering done with the given parameters. The obtained solutions are evaluated on the 4 measures (MDvar, Tvar, ShaP, SPass)
## and then the evolutionary algorithm selects the best solutions for the next generation.
##
## The two evolutionary operators are: crossover and mutation.
## 	crossover - the genome from two parents is used to create an offspring. The values of the parameters are constructed as an
##			average of the values of the two parents. The centroids are chosen randomly between the parents.
##	mutation - randomly mutate the value of a parameter or replace one centroid with a random observation.
##
## Parameters:
##	nind [100] - the number of individuals in each generation;
##	ngmax [100] - maximum number of generations;
##	psurv [0.1] - percentage of survivals between generations. Default is 10% * nind = 0.1 * 100 = 10 individuals
##	pmuts [0.2] - percentage of mutations from the survivals. Default 0.2 * 0.1 * 100 = 2 individuals.
##	maxmut [2] - maximum number of mutants to create (otherwise, no normal children are created).
##	noMutations [2] - number of mutations that a individual suffers when mutating.
##	replacedCentroids [0.1] - percentage of replaced centroids when mutating
##	optimizationCriteria [ [1 2 3 4] ] - a vector with the order number of the criteria to be optimized. Default all. 
##		MDvar = 1, Tvar = 2, ShaP = 3 and SPass = 4. Ex. for optimizing MDvar and ShaP, use [1 3]
##
##	History [ [] ] - History to continue from. History is an internal variable which contains all Configs and all Outputs from all 
##			generations. This option can be used to continue an interrupted execution. History can be either the Octave
##			structure containing the history or the name of the file in which it was saved (see saveHistory).
##	saveHistory ['saved-history.bin.save'] - name of the binary file to save History to. History is saved after each iteration or 
##			after initialization.
##
##	dataset [] - the dataset to be loaded (see temporalKMdeans.m for details).
##	noClusters [5] - number of clusters to construct (see temporalKMdeans.m for details).
##	maxiter [10] - the maximum number of clustering iterations to perform for each individual (see temporalKMdeans.m for details).
##	typeT [5] - type of algorithm to use (see temporalKMdeans.m for details).
##	loadCentroidsFromFile [] - used to load initial centroids for each individual from file (see temporalKMdeans.m for details).
##	identicalInitialCentroids [false] - all individuals start from the same centroids. Children and mutants also (mutations and children
##			do not inherite any centroids. To be used with loadCentroidsFromFile.
##
## WARNING: This script is optimized for multithreading. Please disable the parallelisation in your BLAS implementation by:
##	1) install and select OpenBLAS as your BLAS implementation (details "update-alternatives --display libblas.so.3gf")
##	2) in your shell set: "export OMP_NUM_THREADS=1". Better yet, start octave as: "export OMP_NUM_THREADS=1 ; octave"
##
## Modified: 	07/03/2014 Initial creation
##		13/03/2014 Added the possiblity to optimize only a subset of the criteria
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [History] = genetic_pareto_front( varargin )

	[nind, ngmax, psurv, pmut, noMutations, replacedCentroids, dataset, typeT, noClusters, maxiter, loadCentroidsFromFile, ...
		identicalInitialCentroids, History, saveHistory, optimizationCriteria, maxmut ] = ... 
		process_options (varargin, 'nind', 100, 'ngmax', 100, 'psurv', 0.1, 'pmut', 0.05, 'noMutations', 2, 'replacedCentroids', 0.1, ...
			'dataset', [], 'typeT', 5, 'noClusters', 5, 'maxiter', 10, 'loadCentroidsFromFile', [], 'identicalInitialCentroids', false, ...
			'History', [], 'saveHistory', 'saved-history.bin.save', 'optimizationCriteria', [1 2 3 4], 'maxmut', 2 );

	# sanity checks
	if ( identicalInitialCentroids && isempty(loadCentroidsFromFile) )
		error('Need to specify the initial centroids to be used with for all population!');
	endif

	# print info concerning the optimization criteria
	optimizationCriteria = unique(optimizationCriteria);
	if ( ~isempty(setdiff([1 2 3 4], optimizationCriteria)) )
		printf("Optimizing following criteria: [");
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
	endif

	#variation domain
	#min max for element in the array
	range_vals = [ -1 1 ; ... 	# range for Alpha
		       0 0.0001; ... 	# Beta TODO CPDS1: 0 0.0001; EC: 0 0.05 
		       2 4 ; ... 	# delta TODO CPDS1: 2 4; EC: 0.0001 2 
		       1 1000 ; ...	# lambda1
		       1 1000 ; ... 	# lambda2
		       1 1000 ]; 	# lambda3

	# the number of criteria to optimize: MDvar, Tvar, ShaP and SPass
	nvar= 4 ;
	# calculate the number of individuals that survive and those that mutate
	nsurv= round(psurv * nind);	# nsurv = psurv * nind
	nmut = round(pmut * nsurv);	# nmut = pmut * psurc * nind
	# initial population does full clustering
	initialmaxiter = maxiter; # this makes sure that the initial solutions are not too close to the Pareto front
	
	if ( ~isempty(History) )
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
				warning(sprintf("The given file '%s' is not a valid history file! Starting from scratch."));
				History = [];
			endif
		endif
	endif
	
	# if some centroids were given, load those and give them to each individual in the population.
	initial_centroids = [];
	if ( ~isempty(loadCentroidsFromFile) )
		printf("-> Loading centroids ...");
		[initial_centroids, noClusters] = createCentroidsEval('filetoload', loadCentroidsFromFile);
		centroidfilename = loadCentroidsFromFile ;
		#if at least one character is not readable
		if ( sum(sum(isnan(centroidfilename))) > 0 || rows(centroidfilename) > 1 || min(min( isprint( centroidfilename ) )) == 0 )
			centroidfilename = "<binary-structure>";
		endif
		printf(" done! Loaded %d centroids from file '%s'.\n", noClusters, centroidfilename);
	endif
	
	# check if we need to initialize
	if ( isempty(History) )
		printf("****> Initialization phase...\n");
		# no History, create new one, initialize population and everything
		History = struct() ;
		
		# 1. create the initial population
		# Measures is the variable in which the values for the indicators are kept for a given generation.
		# Each line is an indicator and each column corresponds to an individual.
		Measures = zeros(nvar, nind) ;
	
		# Configs is the structure array which contains the centroids, ctimestamp and clusassign - the output of the run
		# it also contains another structure containing the values of the parameters, which are to be determined: Alpha, Beta, Delta, lamb1, lamb2 and lamb3
		Configs = struct() ;
	
		# and now initialize the first generation at random
		for i=1:nind
			params = struct() ;
			params.Alpha = range_vals(1, 1) + ( range_vals(1, 2) - range_vals(1, 1) ) * rand(1); # init alpha
			params.Beta = range_vals(2, 1) + ( range_vals(2, 2) - range_vals(2, 1) ) * rand(1); # init beta
			params.Delta = range_vals(3, 1) + ( range_vals(3, 2) - range_vals(3, 1) ) * rand(1); # init delta
			params.lamb1 = range_vals(4, 1) + ( range_vals(4, 2) - range_vals(4, 1) ) * rand(1); # init lamb1
			params.lamb2 = range_vals(5, 1) + ( range_vals(5, 2) - range_vals(5, 1) ) * rand(1); # init lamb2
			params.lamb3 = range_vals(6, 1) + ( range_vals(6, 2) - range_vals(6, 1) ) * rand(1); # init lamb3
			# and the fixed part of parameters
			params.maxiter = initialmaxiter ;
			params.dataset = dataset ;
			params.typeT = typeT;
			params.noClusters = noClusters;
		
			Configs(i).params = params ;
			Configs(i).centroids = initial_centroids ;
		endfor

		#=============================== run maxiter iterations for each individual, calculate scores and create the measures matrix
		my_cell = struct2cell(Configs) ;
		printf("--> Evolving individuals using %d clustering iterations ...\n", initialmaxiter);
        ## parallel version
		[Adj, clusassign, centroids_eval, ctimestamp, centroids] = ...
			parcellfun(nproc(), @temporalKMeans_wrapper, my_cell(1,1,:)(:), my_cell(2,1,:)(:), ...
			'UniformOutput', false, 'ErrorHandler', @error_function) ;
        ## sequential version, for debug
#		[Adj, clusassign, centroids_eval, ctimestamp, centroids] = ...
#			cellfun(@temporalKMeans_wrapper, my_cell(1,1,:)(:), my_cell(2,1,:)(:), 'UniformOutput', false) ; # removed , 'ErrorHandler', @error_function
		Outputs = cell2struct( [Adj clusassign centroids ctimestamp], {'Adj', 'clusassign', 'centroids', 'ctimestamp'}, 2) ;
		Outputs = Outputs';
		[MDvar, Tvar, ShaP, SPass] = ...
			parcellfun(nproc(), @evaluateTemporalClustering_wrapper, my_cell(1,1,:)(:), Adj, clusassign, centroids_eval, ctimestamp) ;
	
		Measures = [MDvar' ; Tvar' ; ShaP' ; SPass'] ;
	
		# construction the domination count
		Measures(nvar+1,:) = domination_count( Measures, optimizationCriteria );
		# ordering by domination count
		[foo, index] = sort(Measures(nvar+1, :));
		Measures = Measures(:, index);
		Configs = Configs(index);
		Outputs = Outputs(index);

		ng = 1;
		History(ng).Configs = Configs;
		History(ng).Outputs = Outputs;
		History(ng).Measures = Measures;
	endif

	##======================================== Start iterations
	printf("Saving History to file '%s'... \n", saveHistory);
	save("-binary", saveHistory, "History");
	determine_best_individual('History', History, 'optimizationCriteria', optimizationCriteria, 'generation', ng);
	nonDomin = sum(Measures(nvar+1,:) == 0); #counting non dominated individuals
	while (( nonDomin < nind ) && (ng <= ngmax))
	    printf("\n\n****> Evolving generation %d ..\n", ng+1);
	    # number of nondominated points and survivals from dominated population
	    nsurv = round(psurv * (nind-nonDomin)); # choosing a percentage of dominated survivors            
	    noSurvivals = nonDomin + nsurv;   # see step 3 from article
	    if (noSurvivals == 1)
	    	noSurvivals = 2;
	    endif
	    printf("--> Only the best %d individuals survived ...\n", noSurvivals);
	    
	    # kill less evolved creatures
	    to_kill = [noSurvivals + 1:columns(Measures)];
	    Measures(:, to_kill) = [];
	    Configs(:, to_kill) = [];
	    Outputs(:, to_kill) = [];
	    
	    # for survivals, put the outputed centroids in the Config
	    for i = 1:noSurvivals
	    	if (identicalInitialCentroids)
	    		Configs(i).centroids = initial_centroids ; #same initial centroids
	    	else
		    	Configs(i).centroids = Outputs(i).centroids;
		endif
	    	Configs(i).params.maxiter = maxiter;
	    	# mark survivors in History
	    	History(ng).Outputs(i).survived = true ;
	    endfor
	    
	    # mark extinct in History
	    for i = noSurvivals+1:length(History(ng).Outputs)
	    	History(ng).Outputs(i).survived = false ;
	    endfor
	    
	    # the index of the current child to construct
	    adenf = noSurvivals + 1;

	    ###### creation of mutants
	    numbermutants = min( [floor(pmut*noSurvivals), maxmut] );
	    printf("--> Creating %d mutants ...\n", numbermutants);
	    # first create the centroid pool to choose from if needed
	    centroidPool = [];
	    for i = 1:length(History(ng).Outputs)
	    	centroidPool = [centroidPool ; History(ng).Outputs(i).centroids];
	    endfor
	    for k = 1:numbermutants
	    
	    	# chose a random survivor to mutate
	    	index = round( noSurvivals * rand(1) );
	    	if (index < 1) 
	    		index = 1;
	    	endif
	    	params = Configs(index).params ;
	    	centroids = Outputs(index).centroids;
	    	
	    	for currentMut = 1:noMutations
	    		# generate mutation number
	    		mutationPos = round( 7 * rand(1) );
	    		if (mutationPos < 1) 
	    			mutationPos = 1;
	    		endif
	    		if ((mutationPos == 7) && identicalInitialCentroids)
		    		mutationPos = 6;
		    	endif
	    		
	    		switch mutationPos
	    			case {1}
					params.Alpha = range_vals(1, 1) + ( range_vals(1, 2) - range_vals(1, 1) ) * rand(1); # mutate alpha
				case {2}
					params.Beta = range_vals(2, 1) + ( range_vals(2, 2) - range_vals(2, 1) ) * rand(1); # init beta
				case {3}
					params.Delta = range_vals(3, 1) + ( range_vals(3, 2) - range_vals(3, 1) ) * rand(1); # init delta
				case {4}
					params.lamb1 = range_vals(4, 1) + ( range_vals(4, 2) - range_vals(4, 1) ) * rand(1); # init lamb1
				case {5}
					params.lamb2 = range_vals(5, 1) + ( range_vals(5, 2) - range_vals(5, 1) ) * rand(1); # init lamb2
				case {6}
					params.lamb3 = range_vals(6, 1) + ( range_vals(6, 2) - range_vals(6, 1) ) * rand(1); # init lamb3
				case {7}
					for replacedCentroid = 1:round(replacedCentroids*noClusters)
						# choose randomly a centroid to change
						centroidPos = round( noClusters * rand(1) );
						if (centroidPos < 1) 
							centroidPos = 1;
						endif
						
						# choose randomly an element to replace it
						centroidReplace = round( rows(centroidPool) * rand(1));
						if (centroidReplace < 1) 
							centroidReplace = 1;
						endif

						# replace centroids						
						centroids(centroidPos, :) = centroidPool(centroidReplace, :);
					endfor
					
			endswitch
		endfor
		# and the fixed part of parameters
		params.maxiter = maxiter ;
		params.dataset = dataset ;
		params.typeT = typeT;
		params.noClusters = noClusters;
	    
    	# add the curent mutant to the population
		Configs(adenf).params = params ;
		if (identicalInitialCentroids)
	    		Configs(adenf).centroids = initial_centroids ; #same initial centroids
	    	else
			Configs(adenf).centroids = centroids ;
		endif
		adenf = adenf + 1;
	    endfor
	    
	    ##### creating normal offstrings
	    printf("--> Creating %d normal children ...\n", nind - adenf + 1);
	    while (adenf<=nind)
	    	### a) determine the two parents of the current child
	#         p1=min(1+floor(rand(1)*nonDomin),nonDomin);                               # first parent is a non-dominated parent
	#         p2=min(nonDomin + 1 + floor(rand(1)*nsurv), nonDomin + nsurv);             # second parent is a dominated survivor
		p1 = min(1+floor(rand(1)*noSurvivals),noSurvivals);                               
		p2 = min(1+floor(rand(1)*noSurvivals),noSurvivals);            
		while (p2==p1)                                        
		   # p2=min(nonDomin + 1 + floor(rand(1)*nsurv), nonDomin + nsurv);             # second parent is a dominated survivor
		   p2 = min(1+floor(rand(1)*noSurvivals),noSurvivals);
		endwhile
		
		### b) construct the child's parameters, as a combination of the parameters of the parents
		p1val = struct2cell(Configs(p1).params);
		p2val = struct2cell(Configs(p2).params);
		e = cell(1, rows(range_vals));
		for i = 1:rows(range_vals)
		    alpha = -0.2 + 1.2 * rand(1);      # take a portion from one parent and another franction from another parent
		    e{i} = alpha * p1val{i} + (1-alpha) * p2val{i};
		    if ( e{i} < range_vals(i, 1) )     # verify that I'm still in the limits of my simulation.
		        e{i} = range_vals(i, 1);
		    endif
		    if ( e{i} > range_vals(i, 2) )
		        e{i} = range_vals(i, 2);
		    endif
		endfor
		
		params = cell2struct(e, {'Alpha', 'Beta', 'Delta', 'lamb1', 'lamb2', 'lamb3'}, 2);
		# and the fixed part of parameters
		params.maxiter = maxiter ;
		params.dataset = dataset ;
		params.typeT = typeT;
		params.noClusters = noClusters;
		
		### c) construct the centroids of the child, by randomly selecting them from the centroids of the parents
		p1val = Outputs(p1).centroids;
		p2val = Outputs(p2).centroids;
		centroids = zeros(size(p1val));
		for i=1:rows(centroids)
			chosenparent = round(rand(1));
			if (chosenparent == 0)
				# take centroid from first parent
				choosefrom = p1val;
			else
				# take it from second parent
				choosefrom = p2val;
			endif
			
			if (rows(choosefrom) >= i)
				choosepos = i;
			else
				# means the parent has fewer centroids than i (due to empty centroids at a given moment)
				# just take a random centroid
				choosepos = min(1+floor(rand(1)*rows(choosefrom)), rows(choosefrom)); 
			endif
			
			centroids(i, :) = choosefrom(choosepos, :);
		endfor

		# add the curent child to the population
		Configs(adenf).params = params ;
		if (identicalInitialCentroids)
	    		Configs(adenf).centroids = initial_centroids ; #same initial centroids
	    	else
			Configs(adenf).centroids = centroids ;
		endif
	    	adenf = adenf + 1;
	    endwhile

	    # classement de la nouvelle population
	    #=============================== run maxiter iterations for each individual, calculate scores and create the measures matrix
	    my_cell = struct2cell(Configs) ;
	    # if the initial centroids are the same, then it is useless to do again the survivors
	    if (identicalInitialCentroids)
	    	work_to_do = [noSurvivals+1:length(Configs)]';
	    	Measures(nvar+1, :) = [];
	    else
	    	work_to_do = [1:length(Configs)]';
	    	Outputs = [];
	    	Measures = [];
	    endif
	    
	    printf("--> Evolving individuals using %d clustering iterations ...\n", maxiter);
#    	    [Adj, clusassign, centroids_eval, ctimestamp, centroids] = ...
#    	    	cellfun(@temporalKMeans_wrapper, my_cell(1,1,:)(work_to_do), my_cell(2,1,:)(work_to_do), 'UniformOutput', false); #, 'ErrorHandler', @error_function) ;
	    [Adj, clusassign, centroids_eval, ctimestamp, centroids] = ...
			parcellfun(nproc(), @temporalKMeans_wrapper, my_cell(1,1,:)(work_to_do), my_cell(2,1,:)(work_to_do), ...
			'UniformOutput', false, 'ErrorHandler', @error_function) ;
	    Outputs_new = cell2struct( [Adj clusassign centroids ctimestamp], {'Adj', 'clusassign', 'centroids', 'ctimestamp'}, 2) ;
	    Outputs = [Outputs Outputs_new'];
	    [MDvar, Tvar, ShaP, SPass] = parcellfun(nproc(), @evaluateTemporalClustering_wrapper, my_cell(1,1,:)(work_to_do), Adj, clusassign, centroids_eval, ctimestamp) ;
	    Measures_new = [MDvar' ; Tvar' ; ShaP' ; SPass'] ;
	    Measures = [Measures Measures_new] ;
	    
	    # construction the domination count
	    Measures(nvar+1,:) = domination_count( Measures, optimizationCriteria );
	    # ordering by domination count
	    [foo, index] = sort(Measures(nvar+1, :));
	    Measures = Measures(:, index);
	    Configs = Configs(index);
	    Outputs = Outputs(index);
	   
	    # a new generation
	    ng = ng + 1 ;
	    
	    History(ng).Configs = Configs;
	    History(ng).Outputs = Outputs;
	    History(ng).Measures = Measures;
	    nonDomin = sum(Measures(nvar+1,:) == 0); #counting non dominated survivors
	    printf("Saving History to file '%s'... \n", saveHistory);
	    save("-binary", saveHistory, "History");
	    determine_best_individual('History', History, 'optimizationCriteria', optimizationCriteria, 'generation', ng);
	endwhile # of generation loop

	# find out best point (closest to the origin) from the last generation
	bestIndividual = determine_best_individual('History', History, 'optimizationCriteria', optimizationCriteria, 'generation', length(History));
	
	# generate the 2d cuts of the pareto front
	genetic_graphics(History, bestIndividual);
endfunction

# error function
function y = error_function (s, x)
	error(sprintf("Error in execution of task %d, id %s, message: \n", s.index, s.identifier, s.message)); 
endfunction
