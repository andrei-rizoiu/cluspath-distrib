## Usage: [centroids, newnoclusters] = createCentroidsEval('dataset', dataset, 'noclusters', noclusters, 'filetosave', filetosave, 'filetoload', filetoload )
##
## Function used to initialise the centroids at random. It guaranties that the centroids are different one from the other
## in order to avoid empty clusters. ALSO, might rewrite the number of clusters if there are not enought candidates
##
## Parameters:
##	- dataset - matrix of data individuals from which to choose initial centroids;
##	- noClusters - number of required centroids;
## 	- filenameToSave (optional) - filename to which to dump centroid selection;
##	- filetoload (optional) - in this case, centroids are not constructed at random from the original data, but they are loaded from the file given as parameter.
##
## Examples:
##	[centroids, newnoclusters] = createCentroidsEval('dataset', data, 'noclusters', 5) - choses 5 centroids from the initial dataset data
##	[centroids, newnoclusters] = createCentroidsEval('dataset', data, 'noclusters', 5, 'filetosave', 'initial-centroids.csv') 
##						- same as before and saves to file initial-centroids.csv
##	[centroids, newnoclusters] = createCentroidsEval('filetoload', 'initial-centroids.csv') - loads the centroids from the file 'initial-centroids.csv'
##
## Modified: 	29/10/2013 - adding saving to/loading from file option, created help.

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [centroids, newnoclusters] = createCentroidsEval(varargin)

[dataset, noclusters, filetosave, filetoload] = process_options (varargin, 'dataset', [], 'noclusters', [], 'filetosave', [], 'filetoload', [] );

## check if needed to load from file
if ( ~isempty(filetoload) )
	# load centroids from file
	[foo1, foo2, centroids] = loadCsvFileMissing(filetoload, false);
	# determine how many they are
	newnoclusters = size(centroids, 1);
	
	return;
endif

if ( isempty(dataset) || isempty(noclusters))
	error('Not enough parameters to work with! See help for a list of parameters.');
endif

#how many candidates we want?
maxfactor=100 ; # more than clusters
maxtries = 20 ; # how many times try to complete the centroids before giving up

% now, lets extract our future centroids
indecs = fix (rand(1, noclusters) * rows(dataset));
#dealing with cases of 0
indecs( indecs == 0) = 1;

#now lets deal with doubles
indecs=unique(indecs);
[foo, position] = eliminateDuplicates( dataset( indecs , : ) );
indecs = indecs(position);
tries = 0;
while ( columns(indecs) ~= noclusters)
	tries += 1;
	#we have at least a double

	#generating some more
	newindecs = fix (rand(1, noclusters) * rows(dataset));
	newindecs( newindecs == 0) = 1;

	# adding to the old selection
	indecs = [ indecs, newindecs];
	indecs=unique(indecs);
	[foo, position] = eliminateDuplicates( dataset( indecs , : ) );
	indecs = indecs(position);
	
	#if there are too many, cut them
	if ( columns(indecs) > noclusters)
		indecs = indecs(1, 1:noclusters);
	endif

	# if tried in vain for many times give up (protection from infinite loops)
	if ( tries > maxtries)
		break;
	endif
endwhile

# do the last try
if ( columns(indecs) < noclusters )
	printf("Last resort solution...\n");
	# we are desperate. Put all candidates there
	newindecs = 1:rows(dataset);

	# adding to the old selection
	indecs = [ indecs, newindecs];
	indecs=unique(indecs);
	[foo, position] = eliminateDuplicates( dataset( indecs , : ) );
	indecs = indecs(position);

	#if there are too many, cut them
	if ( columns(indecs) > noclusters)
		indecs = indecs(1, 1:noclusters);
	endif
endif

# still not having them? means there are less candidates then clusters.
# modifying the cluster no

centroids = dataset( indecs , : );
newnoclusters = columns(indecs);

if ( newnoclusters ~= noclusters)
	printf("Not enough candidates... modifying no clusters to %d (from %d)\n", newnoclusters, noclusters );
endif

# if demanded to save to file.
if ( ~isempty(filetosave) )
	# load centroids from file
	saveCsvFile([], [], centroids, filetosave);
endif

endfunction
