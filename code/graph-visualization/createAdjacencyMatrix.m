## Usage: [adj, ctimestamp, arcs, entities, startsat] = createAdjacencyMatrix(files, dtimestamp, clusassign)
##
## Published in the ICTAI 2012 extension in the journal IJAIT 2013.
##
## Based on the entity name, dtimestamp and cluster assignement for each observation, this function
## calculates the static posteriori adjacency matrix between clusassign. It uses a formula based on the time difference and
## the inter_PHI function calculated in the interPhi.m file
##
## Parameters: 	'filename' or ('files', 'dtimestamp', 'clusassign') the data needed
##		'filtcycle' (default false) controls wheather cycles 1 -> 2 -> 1 are moved
##		'omega1' and 'omega2' - are the weights of the temporal distance and the intersection distance between clusters in
##					the adjacency matrix
##		'useDiameters' (defaul true) - if diameters are used (see the "temporalKMeans" function), then we use them to 
##					normalize the temporal distance component. If false, then normalize using the largest value.
##		'ctimestamp' (default []) - the timestamps of centroids. If not timestamp is given (default []), then recompute
##					timestamps by averaging the timestamps of the observations in clusters.
##
## It outputs several values:
##	- adj - the adjacency matrix, where adj(i, j) is the adjacency score between clusters i and j, as proposed in article.
##	- counts - the counts matrix, counts(i,j) is the number of entities passing from cluster i to cluster j
##	- ctimestamp - recalculates the timestamps of centroids;
##	- arcs - arcs{i, j} is a list of entities that pass from state i to state j
##	- entities - are the names of entities, in ascending order;
##	- starts - a vector of the same dimension as the list of entities, so that startsat(i) show the first
##		cluster for entity i.
##
## Example run:
##	[adj, ctimestamp, arcs, countr, startsat] = createAdjacencyMatrix('filename', 'output-tempered.csv', 'filtcycle', false, 'omega1', 1, 'omega2', 4);
##
## Modified:	16/11/2012 Revisited, added this help
##		29/09/2013 Added receiving and passing the ctimestamp

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [adj, counts, arcs, entities, startsat] = createAdjacencyMatrix(varargin)

[filename, files, dtimestamp, clusassign, filtcycle, omega1, omega2, useDiameters, ctimestamp] = ...
    process_options(varargin, 'filename', [], 'files', [], 'dtimestamp', [], 'clusassign', [], 'filtcycle', false, ...
	'omega1', 1, 'omega2', 1, 'useDiameters', true, 'ctimestamp', [] );

	if ( ~isempty(filename) )
		[files, tags, vals] = loadCsvFileMissing(filename, true);
		dtimestamp = vals(:, 1);
		vals(:, 1) = [];
		clusassign = vals(:, end);
		vals(:, end) = [];
	endif
	
	if ( isempty(files) || isempty(dtimestamp) || isempty(clusassign) )
		error("I don't have all the data I need to work!");
	endif

	######################################## ENDED INITIALIZATION ######################

	# first calculate the inter_PHI function value
	# also get the ctimestamp for the temporal difference part
	[interphi, counts, ctimestamp, arcs, entities, startsat] = ...
		interPhi('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', filtcycle, 'ctimestamp', ctimestamp);

	# second calculate the temporal distance between clusters
	# this distance is dt(i,j) = \mu_i^t - \mu_j^t
	mat1 = repmat(ctimestamp, length(ctimestamp), 1); # from which we substract
	mat2 = repmat(ctimestamp', 1, length(ctimestamp)); # what we substract
	dt = abs(mat1 - mat2);
	if (useDiameters)  #calculate the normalization factor, which can be either the temporal diameter or the max value
		normal = calculateDiameter(dtimestamp);
	else
		normal = max(max(dt));
	endif
	dt = dt / normal;

	# get the value of the final measure
	adj = ( omega1 * (1 - dt) + omega2 * ( 1 - interphi) ) / (omega1 + omega2);
endfunction
