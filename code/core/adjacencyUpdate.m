## Usage: [Adj] = adjacencyUpdate('files', [], 'dtimestamp', [], 'ctimestamp', [], 'clusassign', [], 'filtcycle', false, 'lamb2', -1, 'lamb3', -1, 'deltaTmax', -1, 'typeT', 3, 'outputfid', 1, 'missing', false, 'deltaXmax', -1, 'centroids', [], 'Alpha', 0);
##
## Updates the adjacency matrix. It is called in the 4th phase (adjacency matrix update) of the temporalKMeans algorithm.
## It is part of the proposal of detecting typical evolution while infering the adjacency matrix in the clustering process. Only applicable for typeT = 3 (see below)
##
## Returns:
##	'Adj' - the recalculated adjacency matrix. OBS: it returns a null matrix if the typeT is less than 3.
##
## Parameters:
##	'files' - files{i} is the individual to which observation di is associated;
##	'dtimestamp' - the temporal description of observations;
##	'ctimestamp' (default []) - the temporal coordinates of clusters;
##	'clusassign' - the cluster assignement of observations at previous iteration;
##	'filtcycle' (default false) controls wheather cycles 1 -> 2 -> 1 are moved;
##	'typeT' (optional, default 3) - selects type of algorithm to use. Description in "temporalKMeans.m";
##	'lamb2', 'lamb3' - when typeT=3, these are the weights of the three components of the objective function
##				(temporal distance of clusters, intersection of clusters in terms of entitites);
##	'deltaXmax' (default -1) - the spatial diameter. -1 means no diameter, the max in the matrix is taken
##	'deltaTmax' (optional, default -1) - the temporal diameter. -1 means no diameter, the max in the matrix is taken
##	'outputfid' (default 1) - the FID to which to print debug messages. 1 is the normal output, 2 is the error, 0 is silent (no debug) 
##	'missing' (default false) - does your data contain missing values? If missing data is detected, that it is set automatically to true.
##	'centroids' (default []) - the centroids at the current iteration
##	'Alpha' (default 0) - the temporal cluster cohesion penalisation factor weight. -1 only temp component, 0 equal importance, 1 only multidimensional component.  
##	'Beta' (default 1) - the temporal entity cohesion penalisation factor weight. Set to zero for no cohesion information (only for typeT = 5)
##	'Delta' (default 3) - the width of the contiguous segmentation function. Set Delta = 3 in order to have around 0.4Beta at a DeltaT=4 (only for typeT = 5)
##
## Modified:	29/09/2013 Created first version.
##		01/11/2013 Added support for typeT = 4
##		27/11/2013 Added support for typeT = 5
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [Adj] = adjacencyUpdate(varargin)

	[files, dtimestamp, ctimestamp, clusassign, filtcycle, lamb2, lamb3, deltaTmax, typeT, outputfid, missing, deltaXmax, centroids, Alpha, Beta, Delta] = ...
		process_options(varargin, 'files', [], 'dtimestamp', [], 'ctimestamp', [], 'clusassign', [], 'filtcycle', false, ...
		'lamb2', -1, 'lamb3', -1, 'deltaTmax', -1, 'typeT', 3, 'outputfid', 1, 'missing', false, 'deltaXmax', -1, ...
		'centroids', [], 'Alpha', 0, 'Beta', 1, 'Delta', 3);
		
	if ( isempty(files) || isempty(dtimestamp) || isempty(clusassign) || isempty(ctimestamp) || (lamb2 == -1) || (lamb3 == -1))
		error("I don't have all the data I need to work!");
	endif
	
	# this function does not work with typeT < 3
	if (typeT < 3)
		Adj = zeros(length(ctimestamp));
		return;
	endif
	
	if ( (sum([4 5] == typeT) > 0) && isempty(centroids))
		error('typeT = 4 or 5, adjacency matrix update requires the centroids parameter!');
	endif
	
	######### Calculate the K matrix
	ctdiff = [];
	# this term is issued from the temporal penalisation from the contiguity semi-supervised penalisation. Only for typeT = 5
	temp_pen = [];
	switch typeT
		case 3
			# the first term T1 is based on the temporal similarity of centroids
			tm1 = repmat( ctimestamp, length(ctimestamp), 1);
			tm2 = repmat( ctimestamp', 1, length(ctimestamp));
			ctdiff = (tm1 - tm2) .^ 2;
			if (deltaTmax > -1)
				divv = deltaTmax .^ 2;
			else
				divv = max(max(ctdiff));
			endif
			ctdiff = ctdiff / divv;
		
		case 4
			# we use the temporal-aware similarity between centroids
			[ctdiff, foo] = distanceTemporal('Amultidim', centroids', 'Atimestamp', ctimestamp, 'Bmultidim', centroids', 'Btimestamp', ctimestamp, ...
				'Alpha', Alpha, 'typeT', typeT, 'outputfid', outputfid, 'missing', missing, 'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, ...
				'outputfid', 0 );
				
		case 5
			# the term containing the temporal-aware dissimilarity between centroids is the same as for typeT = 4
			# we use the temporal-aware similarity between centroids
			[ctdiff, foo] = distanceTemporal('Amultidim', centroids', 'Atimestamp', ctimestamp, 'Bmultidim', centroids', 'Btimestamp', ctimestamp, ...
				'Alpha', Alpha, 'typeT', typeT, 'outputfid', outputfid, 'missing', missing, 'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, ...
				'outputfid', 0 );
				
			# calculate the temporal penalty factor
			temp_pen = zeros(size(ctdiff));
			noclusters = rows(temp_pen);
			
			# for each pair of clusters (r, s), we need to determine all tranzitions and associated temporal cost
			for r=1:noclusters
				for s=1:noclusters
					if (r ~= s)
						members_r = find(clusassign == r); # index of observations in cluster r
						members_s = find(clusassign == s); # index of observations in cluster s
					
						# determine if the same entity
						# ent_r(i, :) have all the same value, the entity to which belongs x_i
						ent_r = repmat(files(members_r)', 1, length(members_s)); 
						# ent_s(:, j) have all the same value, the entity to which belongs x_j
						ent_s = repmat(files(members_s), length(members_r), 1);  
					
						select = strcmp(ent_r, ent_s); # 1 for same entity, 0 for different entities
						# apparently, strcmp returns a boolean matrix and not an integer matrix. Need to transform it before applying 
						# the NaN setting
						select = double(select) ;
						select(select == 0) = NaN; # so that, when multiplying, it makes NaN for any operation.
					
						# determine precedence and temporal difference. Remember: dtimstamp is a column vector
						tdiff_r = repmat(dtimestamp(members_r), 1, length(members_s)); 
						tdiff_s = repmat(dtimestamp(members_s)', length(members_r), 1); 
						
						tdiff = tdiff_r - tdiff_s;
						tdiff(tdiff >= 0) = NaN;   # if tdiff >=0 means that x_i^t >= x_j^t, with x_i in r and x_j in s
					
						# apply selection
						tdiff = tdiff .* select;
						
						# calculate the temporal penalty
						# we force the normal distribution-inspired penalty function
						tdiff = Beta * exp( -0.5 .* ((tdiff / Delta) .^ 2) );
					
						# sum everything, except NaN's into the penalty term of (r, s)
						temp_pen(r, s) = sum(sumskipnan(tdiff)); 
					else
						temp_pen(r, s) = 0;
					endif
				endfor
			endfor
				
		otherwise
			# can only be higher that 4, since for typeT<3 it returns an empty matrix at the verification phase
			error(sprintf("I don't recognize 'typeT=%d' Did you give me a wrong typeT algorithm to work with???", typeT));
	endswitch
	
	T1 = lamb2 * ctdiff;
	
	# calculate second term T2, the inter_phi intersection between clusters 
	[interphi] = interPhi('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', filtcycle, 'ctimestamp', ctimestamp);
	interphi = interphi .^ 2;
	T2 = lamb3 * interphi;
	
	# finally, the temporal penalisation term. If temp_pen is empty, then we are in typeT < 5 and set it to zero
	if (~isempty(temp_pen))
		T3 = temp_pen;
	else
		T3 = zeros(size(T1));
	endif
	
	K = T1 + T2 - T3;
	
	######## Calculate the fixed term sum sum 1 / Kpq
	# put NaNs on the main diagonal of K, since we do not want to calculate clusters with themselves
	nullify = diag(repmat(NaN, 1, rows(K))) + 1; # a matrix of ones, having NaN on main diagonal
	K = K .* nullify;
	TMP = 1 ./ K;
	
	FT = sumskipnan(sumskipnan(TMP));
	
	######## Finally construct the matrix
	Adj = K .* FT;
	Adj = 1 ./ Adj;
	Adj( isnan(Adj) ) = 0;	#replace the NaN's on the main diag with zeros
	
endfunction
