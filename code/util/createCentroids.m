function [centroids, newnoclusters] = createCentroids(fn, noclusters)

# function used to initialise the centroids at random.
# it guaranties that the centroids are different one from the other
# in order to avoid empty clusters
# ALSO, might rewrite the number of clusters if there are not enought candidates

#how many candidates we want?
maxfactor = 5 ; # more than clusters
minfiles = 10 ; # candidates from minimum minfiles
maxtries = 20 ; # how many times try to complete the centroids before giving up

if ( minfiles > columns(fn) )
	minfiles = columns(fn);
endif

#load the minfiles
matrix = [];
for i=1:minfiles
	matrix = [matrix ; loadFile( fn{i} ) ];
endfor

#if there are not enough candidates, continue until ratio
fileno = minfiles + 1;
while ( rows(matrix) < maxfactor * noclusters)
	if (fileno <= columns(fn) )
		matrix = [matrix ; loadFile( fn{fileno} ) ];
		fileno += 1;
	else
		# no more files to load, so we stop anyways
		break;
	endif
endwhile

printf("Selected candidates from %d files to create %d centroids...\n", fileno-1, noclusters );

% now, lets extract our future centroids
indecs = fix (rand(1, noclusters) * rows(matrix));
#dealing with cases of 0
indecs( indecs == 0) = 1;

#now lets deal with doubles
indecs=unique(indecs);
[foo, position] = eliminateDuplicates( matrix( indecs , 6:end ) );
indecs = indecs(position);
tries = 0;
while ( columns(indecs) ~= noclusters)
	tries += 1;
	#we have at least a double

	#generating some more
	newindecs = fix (rand(1, noclusters) * rows(matrix));
	newindecs( newindecs == 0) = 1;

	# adding to the old selection
	indecs = [ indecs, newindecs];
	indecs=unique(indecs);
	[foo, position] = eliminateDuplicates( matrix( indecs , 6:end ) );
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
	newindecs = 1:rows(matrix);

	# adding to the old selection
	indecs = [ indecs, newindecs];
	indecs=unique(indecs);
	[foo, position] = eliminateDuplicates( matrix( indecs , 6:end ) );
	indecs = indecs(position);

	#if there are too many, cut them
	if ( columns(indecs) > noclusters)
		indecs = indecs(1, 1:noclusters);
	endif
endif

# still not having them? means there are less candidates then clusters.
# modifying the cluster no

centroids = matrix( indecs , 6:end );
newnoclusters = columns(indecs);

if ( newnoclusters ~= noclusters)
	printf("Not enough candidates... modifying no clusters to %d (from %d)\n", newnoclusters, noclusters );
endif

end
