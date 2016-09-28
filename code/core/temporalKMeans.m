## Usage: 
##	************* for typeT = 0 ************************
## [clusassign, dtimestamp, centroids, ctimestamp] = temporalKMeans('file', 'temporal-average.csv', 'noClusters', 5, 'Alpha', 1.5, 'Beta', 2, 'Delta', 3, 'replaceMissing', false, 'normalize', true, 'stardardize', true, 'maxiter', 100, 'typeT', 0, 'saveResults', true);
##
##	************* for typeT = 1 ************************
## [clusassign, dtimestamp, centroids, ctimestamp] = temporalKMeans('file', 'temporal-average.csv', 'noClusters', 10, 'Alpha', 1.5, 'Beta', 0.003, 'Delta', 3, 'replaceMissing', false, 'normalize', true, 'stardardize', true, 'maxiter', 100, 'typeT', 1, 'saveResults', true);
##	************* for typeT = 2 ************************
## [clusassign, dtimestamp, centroids, ctimestamp] = temporalKMeans('file', 'temporal-average.csv', 'noClusters', 10, 'Alpha', 1.5, 'Beta', 3, 'Delta', 5, 'replaceMissing', false, 'normalize', true, 'stardardize', true, 'maxiter', 100, 'typeT', 2, 'saveResults', true);
##
## An implementation of KMeans used in regrouping the documents. Centroids are initialized at random.
##	'file' - the data file to be loaded with "loadTemporalCsv" function;
##	'noClusters' (optional, default 5) - no of groups in whith to partition;
##	'Alpha' (optional default 0 for typeT=0, 0 for typeT=1) - the temporal cluster cohesion penalisation factor weight. 
##			-1 only temp component, 0 equal importance, 1 only multidimensional component (Euclidean distance). Proposed in IJAIT 13
##	'Beta' (optional default 0) - the temporal individual cohesion penalisation factor weight. Set to zero for no cohesion information.
##	'Delta' (optional, default 3) - the width of the contiguous segmentation function;
##	'replaceMissing' (optional, default false) - replace missing values with mean or just ignore them when calculating the measure;
##	'normalize' (optional, default true) - should we normalize the dataset before performing calculations?
##	'standardize' (optional, default true) - should we standardize instead of normalizing?
##	'maxiter' (optional, default 50) - maximum no of iterations;
##	'threshold' (optional, default 0.001) - the absolute difference in objective function when to stop
##	'typeT' (optional, default 1) - selects type of algorithm to use. Description below.
##	'saveResults' (optional, default false) - if true, the results (the dataset with the cluster assigned, the centroids and various graphics) are saved to disk. 
##						If false, nothing is saved, only the return values are returned.
##	'outputfid' (optional, default 1) - the FID to which to print debug messages. 1 is the normal output, 2 is the error, 0 is silent (no debug)
##	'useDiameters' (optional, default true) - if set, use temporal and spatial diameters when calculating distances, otherwish use the maximum in the matrix 
##						(NOT RECOMMENDED TO CHANGE - just leave default)
##	'lamb1', 'lamb2', 'lamb3' (all default 1) - when typeT=3 (see after), these are the weights of the three components of the objective function
##						(assignement of observations to clusters, temporal distance of clusters, intersection of clusters in terms of entitites).
##						All lamb are normalized so that their sum is 1 (default, they all get to be 0.33) (NOT ANY MORE)
##	'loadCentroidsFromFile' (default []) - load the initial centroids from file instead of initializing them at random from dataset. Usefull for comparing different
##						algorithms on the same initial centroids.
##	'saveInitCentroids' (default true) - when centroids are generated, do I save them to the file 'init-centroids.csv' ?
##
##	'typeT' definitions:
##		- typeT = 0 - (OBSOLETE, NEVER PUBLISHED) uses Euclidian distance for calculating the distance in the multimensional 
##				space. Uses a cannot-link component to give the clusters consistence in the temporal space. Uses the 
##				must-link component to assure contiguous segmentation for entities. (function inspired from normal)
##		- typeT = 1 - (ICTAI 2012) temporal-aware dissimilarity function for multidimensional and temporal cohesion.
##				Same must-link contiguous segmentation function as 0
##		- typeT = 2 - (LIN06 - in Aigaion). Euclidean distance, no temporal cohesion, window based contiguous 
##				segmentation function.
##		- typeT = 3 - (UNPUBLISHED) algorithm infering graph structure inside the clustering (as proposed initially 
##              by Julien Velcin in the 2013.06 work document)
##		- typeT = 4 - (UNPUBLISHED) similar to typeT = 3, but the second term of the objective function is now using 
##              the temporal-aware distance between centroids
##		- typeT = 5 - (current) similar to typeT = 4, but the temporal contiguity penalty function has been modified 
##              to take into account the adjacency matrix
##
## Modified: 	
##      14/11/2012 Put default maxiter parameter to 50 too long execution time.
##		27/09/2013 Improved the help
##		29/09/2013 Added support for typeT = 3
##		01/11/2013 Added support for typeT = 4
##		27/11/2013 Added support for typeT = 5
##      01/12/2014 Fixed a bug when all but one cluster was eliminated in the clustering process
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [Adj, clusassign, dtimestamp, centroids, ctimestamp, centroids_EVOLUTION] = temporalKMeans( varargin )

############################################ Read parameters ################################

[file, noClusters, Alpha, Beta, Delta, replaceMissing, normalize, standardize, maxiter, threshold, typeT, saveResults, outputfid, ... 
useDiameters, lamb1, lamb2, lamb3, loadCentroidsFromFile, saveInitCentroids ] = ... 
	process_options (varargin, 'file', [], 'noClusters', 5, 'Alpha', 0, 'Beta', 0, 'Delta', 3, 'replaceMissing', false, ...
		'normalize', true, 'standardize', true, 'maxiter', 50, 'threshold', 0.001, 'typeT', 3, 'saveResults', false, ... 
		'outputfid', 1, 'useDiameters', true, 'lamb1', 1, 'lamb2', 1, 'lamb3', 1, 'loadCentroidsFromFile', [], 'saveInitCentroids', true );

if ( isempty(file) )
	error('I need the file to work with!');
endif

########################################### Some Variables ##############################################

missing = false;   #do we have missing data in the dataset

# if we don't want another algorithm then typeT = 3, then lamb1 = 1 and the others are 0
if ( typeT < 3 )
	lamb1 = 1;
	lamb2 = 0;
	lamb3 = 0;
endif

########################################### Pretreatment ################################################

fprintf(outputfid, "--> Parameters: Alpha = %.5f, Beta = %.5f, Delta = %.5f, lamb1 = %.5f, lamb2 = %.5f, lamb3 = %.5f\n", ...
                    Alpha, Beta, Delta, lamb1, lamb2, lamb3 );

fprintf(outputfid, "-> Load data file ...");
[files, tags, data] = loadCsvFileMissing ( file, true );
fprintf(outputfid, " done!\n");

# isolate the timestamp from the data
dtimestamp = data(:, 1);
data(:, 1) = [];	#eliminate timestamp information

# if demanded, replace missing values by averages
if ( replaceMissing )
	fprintf(outputfid, "-> Replacing missing values by means ...");
	data = replaceMissingByMean( data );
	fprintf(outputfid, " done!\n");
	missing = false;
else
	if ( sum(sum(isnan(data))) > 0)
		#means we have missing data
		missing = true;
	endif
endif

#make a copy of initial data and normalize if necessary
Dinitial = data;
DtimestampInitial = dtimestamp;
if (normalize)
	fprintf(outputfid, "-> Normalizing data ...");
	data = normalizeMatr(Dinitial, standardize);

	# normalising timestamps
	# this is necessary in order not to put too much weight on the temporal component relative to the descriptive component
	dtimestamp = normalizeMatr(DtimestampInitial, standardize);
	
	# save new stardadized file to disk
	saveCsvFile(files, tags, [dtimestamp data], "input-normalized.csv");

	# done
	fprintf(outputfid, " done!\n");
endif

deltaTmax = -1;
deltaXmax = -1;
if (useDiameters)
	fprintf(outputfid, "-> Calculating diameters ...");
	deltaTmax = calculateDiameter(dtimestamp, missing);
	deltaXmax = calculateDiameter(data, missing);
	fprintf(outputfid, " done!\n");
endif

########################################### Start algo ################################################

#initialize main variables
dataDim = columns(data);
noData = rows(data);

restartCount = 0;
restartMax = 50; # we don't wanna restart too many times either
restart = true;
# keep everything in a loop, so that we can restart
while ( restart && (restartCount < restartMax) )
	# by default only once
	restart = false;
	restartCount = restartCount + 1;
	if (restartCount > 1)
		fprintf(outputfid, "--> Some error occured, try %d/%d...\n\n", restartCount, restartMax);
	endif
	
	% init the centroids randomly
	fprintf(outputfid, "-> Initializing centroids ...");
	# put back timestap so that centroids have timestaps attached
	dataT = [dtimestamp data];
	if ( ~isempty(loadCentroidsFromFile) )
		[centroids, noClusters] = createCentroidsEval('filetoload', loadCentroidsFromFile);
		centroidfilename = loadCentroidsFromFile ;
		#if at least one character is not readable
		if ( sum(sum(isnan(centroidfilename))) > 0 || rows(centroidfilename) > 1 || min(min( isprint( centroidfilename ) )) == 0 )
			centroidfilename = "<binary-structure>";
		endif
		fprintf(outputfid, " done! Loaded %d centroids from file '%s'.\n", noClusters, centroidfilename);
	else
		if (saveInitCentroids)
			[centroids, noClusters] = createCentroidsEval('dataset', dataT, 'noClusters', noClusters, 'filetosave', 'init-centroids.csv');
		else
			[centroids, noClusters] = createCentroidsEval('dataset', dataT, 'noClusters', noClusters);
		endif
		fprintf(outputfid, " done! Constructed %d centroids.\n", noClusters);
	endif
	ctimestamp = centroids(:, 1)';
	centroids(:, 1) = [];	#eliminate timestamp information

	fprintf(outputfid, "-> Starting regrouping iterations ...\n");
	%repeat until the partition doesn't change anymore
	myDiverg = 0.0;
	oldDiverg = inf;
	iter = 0;
	clusassign = [];
	Adj = zeros(noClusters, noClusters);	# initialize the adjacency matrix to zero.

	while ( abs(oldDiverg - myDiverg) > threshold )
		%starting time
		tic

		if ( iter > 0 )
			oldDiverg = myDiverg;
		endif
		iter++;

		result = zeros(noClusters, dataDim);		#used, together with division, to update centroids in M phase
		division = zeros(noClusters, dataDim);
		minim = [];	#minim(i) is the temporal distance between observation i and its cluster clusAssig(i) (including temporal info)

		###### PHASE 1 - calculating the distance between individuals and the centroids

		# if the cluster assignement is empty, we cannot use our temporal information
		# so do an empty round, just to calculate assignements
		if ( isempty(clusassign) && (sum([1 3 4 5] == typeT) > 0) )
			fprintf(outputfid, "--> Temporal Distance Calculation (extra step) ... ");
			[DTemp, D] = distanceTemporal('Amultidim', data', 'Atimestamp', dtimestamp, 'Bmultidim', centroids', ...
				'Btimestamp', ctimestamp, 'files', files, 'clusassign', clusassign, 'Alpha', Alpha, ...
				'Beta', Beta, 'Delta', Delta, 'typeT', typeT, 'outputfid', outputfid, 'missing', missing, ...
				'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, 'Adj', Adj );
			[minim, clusassign] = min (DTemp');
		endif
		fprintf(outputfid, "--> Temporal Distance Calculation ... ");
		[DTemp, D] = distanceTemporal('Amultidim', data', 'Atimestamp', dtimestamp, 'Bmultidim', centroids', ...
				'Btimestamp', ctimestamp, 'files', files, 'clusassign', clusassign, 'Alpha', Alpha, ...
				'Beta', Beta, 'Delta', Delta, 'typeT', typeT, 'outputfid', outputfid, 'missing', missing, ...
				'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, 'Adj', Adj );
		
				
		###### PHASE 2 - calculating the assignement of individuals to clusters
		
		[minim, clusassign] = min (DTemp');

		# calculate square Euclidean distance between observation i and cluster clusassign(i)
		# this is needed in the centroid multidimensional update phase, when diameters are not used.
		fprintf(outputfid, "--> Assign observations to clusters ...\n");
		distObsClust = zeros(1, noData);	# square eucludian distance between observation i and its cluster, clusAssig(i)
		for i=1:noData
			distObsClust(i) = D(i, clusassign(i));
		endfor
		
		###### PHASE 3 - update the centroids
				
		% for the stoppingCriterion
		oldPositions = centroids;
		oldCtimestamp = ctimestamp;

		fprintf(outputfid, "--> Update centroids ...\n");		
		[centroids, ctimestamp, pointsInCluster] = centroidUpdate('data', data, 'dtimestamp', dtimestamp, ...
				'distObsClust', distObsClust, 'clusassign', clusassign, 'typeT', typeT, ...
				'noClusters', noClusters, 'oldCtimestamp', oldCtimestamp, 'Alpha', Alpha, ...
				'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, 'lamb1', lamb1, ...
				'lamb2', lamb2, 'Adj', Adj, 'oldCentroids', centroids );
				
		#eliminate empty clusters
		if ( sum (pointsInCluster == 0) > 0 )
			# there are empty clusters, we should eliminate them
			empt = (pointsInCluster == 0);
			# which ones are they? after the next instruction, empts will contain the index of empty clusters
			empts = find(empt);
			# remove from old cluster records
			oldPositions(empts, :) = [];
			# remove correponding rows and columns from the Adj matrix
			Adj(empts, :) = [];
			Adj(:, empts) = [];
			
			# need to eliminate them from clusassign and degrade all others
			for i=1:length(empts)
				eclus = empts(i);
				clusassign(clusassign == eclus) = [];	#normally this should never happen, it is an empty cluster
				#we eliminated one, so the order of all the others should be decreased by one
				clusassign(clusassign > eclus)--; 
				centroids(eclus, :) = [];
				ctimestamp(eclus) = [];
				noClusters--;
				
				empts(empts>eclus)--;
			endfor
			
			fprintf(outputfid, "-> WE HAVE %d EMPTY CLUSTER(s)!!!!!! Eliminating them!\n", sum(empt));
		endif
		
		###### PHASE 4 - update adjacency matrix
		
		fprintf(outputfid, "--> Update adjacency matrix ...\n");
		[Adj] = adjacencyUpdate('files', files, 'dtimestamp', dtimestamp, 'ctimestamp', ctimestamp, 'clusassign', clusassign, ...
				'typeT', typeT, 'filtcycle', false, 'lamb2', lamb2, 'lamb3', lamb3, 'deltaTmax', deltaTmax, ...
				'outputfid', outputfid, 'missing', missing, 'deltaTmax', deltaTmax, 'centroids', centroids, 'Alpha', Alpha, ...
				'Beta', Beta, 'Delta', Delta );
		
		###### PHASE 5 - calculate objective measure
		# part coresponding to the first term (T1) in the objective function: determined by the assignement of observations to clusters
		foptT1 = sum(minim);
		
		# part coresponding to the second term (T2) in the objective function: stating that clusters should be close in time
		switch typeT
			case 3
				######################### this is for typeT = 3 ###############################
				t1 = repmat( ctimestamp, length(ctimestamp), 1);
				t2 = repmat( ctimestamp', 1, length(ctimestamp));
				ctdiff = (t1 - t2) .^ 2;
				if (useDiameters)
					divv = deltaTmax .^ 2;
				else
					divv = max(max(ctdiff));
				endif
				ctdiff = ctdiff ./ divv;
		
			case 4
				######################### switched to typeT = 4 ################################
				[ctdiff, D] = distanceTemporal('Amultidim', centroids', 'Atimestamp', ctimestamp, 'Bmultidim', centroids', 'Btimestamp', ctimestamp, ...
						'Alpha', Alpha, 'typeT', typeT, 'outputfid', outputfid, 'missing', missing, 'deltaXmax', deltaXmax, 'deltaTmax', ... 
						deltaTmax );
			otherwise
				ctdiff = zeros(size(Adj));
		endswitch
		ctdiff = (Adj .^ 2) .* ctdiff ;
		
		foptT2 = sum(sum(ctdiff));
		
		# part coresponding to the third term (T3) in the objective function: intersection of clusters in terms of entities.
		[interphi] = interPhi('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', false, 'ctimestamp', ctimestamp);
		interphi = (Adj .^ 2) .* (interphi .^ 2);
		foptT3 = sum(sum(interphi));
		
		# and now sum it up
		myDiverg = lamb1*foptT1 + lamb2*foptT2 + lamb3*foptT3;
		fprintf(outputfid, "Means: T1: %.3f, T2: %.3f, T3: %.3f\n", lamb1*foptT1, lamb2*foptT2, lamb3*foptT3);

		iter_time = toc;
		fprintf(outputfid, "Iteration %d, OV %.3f, dif: %.3f - time: %.2f sec\n", iter, myDiverg, oldDiverg - myDiverg, iter_time);
		
		%stoppingCriterion
		pos_diff = sum (sumskipnan( (centroids - oldPositions) .^2 ) );

        # check if maximum number of iterations was already performed
		if ( iter >= maxiter ) 
			break;
		endif
		
		# sanity check - what if there is only one cluster left. Then there is nothing more we can do.
		if (noClusters == 1)
			fprintf(outputfid, "***> Only one cluster left! Ending clustering. Check your parameters and/or initial centroids!\n");
			break;
		endif

	endwhile

endwhile

# the next return value are the centroids which are needed in case we want to continue the evolution
# they are useful for the evolutionary algorithm
centroids_EVOLUTION = [vec(ctimestamp) centroids];

# sort the centroids so that they will be in chronological order
[foo, index] = sort(ctimestamp);
ctimestamp = foo;
centroids = centroids(index, :);
[foo, indexR] = sort(index);
clusassign = indexR(clusassign);
Adj = Adj(index,index);

# restore the timestamps of data and rebuild those of centroids
dtimestamp = DtimestampInitial;
ctimestamp = deNormalizeMatr(ctimestamp', dtimestamp, standardize);
ctimestamp = ctimestamp';

# if needed, reconstruct adjacency matrix as proposed in IJAIT 14
if (typeT < 3)
	Adj = createAdjacencyMatrix('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', false, 'omega1', 0, 'omega2', 1, 'useDiameters', useDiameters, 'ctimestamp', ctimestamp);
endif

#save the results
if (saveResults)
	[files, tags, data] = addClusterInfo(file, clusassign, 'output-clustered.csv');
	# first column of data is the timestamp of observations and last column of data is the cluster assignement
	tcentroids = [vec(ctimestamp) centroids];
	cfiles = {};
	for i = 1:rows(centroids)
		cfiles = [cfiles {sprintf("cluster%d", i)} ];
	endfor

	saveCsvFile(cfiles, tags, tcentroids, "output-centroids.csv");
	
	# construct the graphic representation of the graph
	switch typeT
		case {0, 1, 2}
			# create the graph file only with interphi component (omega=1) and no filtering
			# another option is to leave by default: omega=-1 (equal importance to temporal and interphi) and 
			# filter by Adj(1,1)
			generateGraphFile('files', files, 'tags', tags, 'vals', data, 'outputfile', 'output-graph.dot', ... 
				'edgedescription', false, 'omega', 1, 'ctimestamp', ctimestamp, 'graph_threshold', 0);
		case {3, 4, 5}
			# give the drawing function the adjacency matrix to plot. We need to set the threshold for filtering
			# set it using the k-1 heuristic
			generateGraphFile('files', files, 'tags', tags, 'vals', data, 'outputfile', 'output-graph.dot', ... 
				'edgedescription', false, 'Adj', Adj, 'ctimestamp', ctimestamp ); 
	endswitch
	
	# set off the visibility to speed up rendering
	set(0, 'defaultfigurevisible', 'off');
	descriptiveGraphics('files', files, 'tags', tags, 'vals', data, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'saveToFile', true, 'dumpEpsFiles', true);
	
	# close all created figures
	close all;
endif

endfunction


