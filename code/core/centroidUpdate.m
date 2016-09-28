## Usage: [centroids, ctimestamp, pointsInCluster] = centroidUpdate('data', data, 'dtimestamp', dtimestamp, 'distObsClust', distObsClust, 'clusassign', clusassign, 'typeT', typeT, 'noClusters', noClusters, 'oldCtimestamp', oldCtimestamp, 'deltaXmax', deltaXmax, 'deltaTmax', deltaTmax, 'lamb1', lamb1, 'lamb2', lamb2, 'Adj', Adj );
##
## Updates the centroids multidimensional descriptive component and the temporal component. It is called in the centroid update phase of the temporalKMeans algorithm.
## Algorithm typeT = 1 was proposed in ICTAI 2012 and typeT = 3 is under construction.
##
## Returns:
##	'centroids' - a matrix with the new spatial component for the centroids;
##	'ctimestamp' - the new temporal component for the centroids;
##	'pointsInCluster' - returns how many observations are in each cluster (useful for detecting empty clusters).
##
## Parameters:
##	'data' - the multidimensional description of observations;
##	'dtimestamp' - the temporal description of observations;
##	'clusassign' - the new assignement of observations into clusters, calculated in this iteration, at the previous step;
##	'typeT' (default 0)- the algorithm type (see "temporalKMeans.m" file for details);
##	'noClusters' - the number of clusters;
##	'distObsClust' (default []) - distObsClust(i) is the SQUARE of the Euclidean distance between observation i and its assigned centroid
##	'oldCtimestamp' (default []) - the OLD temporal coordinates of clusters, from the previous iteration. OBS: the spatial coordinates of 
##			clusters are embedded in 'distObsClust'
##	'Alpha' (default 1) - the temporal cluster cohesion penalisation factor weight. 
##			-1 only temp component, 0 equal importance, 1 only multidimensional component (Euclidean distance). Proposed in IJAIT 13
##	'deltaXmax' (default -1) - spatial diameter of the dataset. If -1, use the largest in local collection;
##	'deltaTmax' (default -1) - temporal diameter of the dataset. If -1, use the largest in local collection;
##	'lamb1', 'lamb2' [only for typeT=3] -  these are the weights of the components of the objective function
##			(assignement of observations to clusters, temporal distance of clusters).
##	'Adj' [only for typeT>=3] - the adjacency matrix calculated during the clustering.
##	'oldCentroids' (only for typeT=4) - the old values for the centroids
##
## Modified:	29/09/2013 Added this help, adapted for typeT = 3 (graph structure inside clustering);
##		01/11/2013 Added support for typeT = 4

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [centroids, ctimestamp, pointsInCluster] = centroidUpdate(varargin)

	[data, dtimestamp, clusassign, typeT, noClusters, distObsClust, oldCtimestamp, Alpha, deltaXmax, deltaTmax, lamb1, lamb2, Adj, oldCentroids] = ...
		process_options(varargin, 'data', [], 'dtimestamp', [], 'clusassign', [], 'typeT', 0, 'noClusters', [], ...
		'distObsClust', [], 'oldCtimestamp', [], 'Alpha', 1, 'deltaXmax', -1, 'deltaTmax', -1, 'lamb1', -1, ...
		'lamb2', -1, 'Adj', [], 'oldCentroids', [] );

	if (isempty(data) || isempty(dtimestamp) || isempty(clusassign) || isempty(noClusters) )
		error('I need at least my data to work with');
	endif
	
	if (sum([3 4] == typeT) > 0)
		if ( (lamb1 == -1) || (lamb2 == -1) || (isempty(Adj)) )
			error('For this algorithm, we need to set extra parameters: lamb1, lamb2 and Adj. (check help)');
		endif
	endif
	
	if (sum([4] == typeT) > 0)
		if ( isempty(oldCentroids) )
			error('For this algorithm, we need to set extra parameters: oldCentroids (check help)');
		endif
	endif

	############################ END sanity check, start algorithm #####################################

	pointsInCluster = zeros(1, noClusters);		#for each cluster, how many individ are in it
	centroids = zeros(noClusters, columns(data) );
	ctimestamp = zeros(1, noClusters);
	# calculate gammaD and gammaT
	gammaD = 1 + Alpha ; 
	if (gammaD > 1)
		gammaD = 1;
	endif

	gammaT = 1 - Alpha ; 
	if (gammaT > 1)
		gammaT = 1;
	endif

	switch typeT
	
		############# perform a simple Euclidean update, like the classical KMeans ###########################	
		case {0, 2}
		for i=1:noClusters
			a =  (clusassign == i); 	#which individ belong in this cluster
			pointsInCluster(i) = sum(a);	#how many individ are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				[rez, div] = sumskipnan( data( a, :) );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs

				% the new centroids will be the average of all the elements in their group
				centroids( i , : ) = rez ./ div;

				[ctimestamp(i), cnt] = sumskipnan( dtimestamp(a) );
				if (cnt>0)
					ctimestamp(i) = ctimestamp(i) / cnt;
				endif
			endif
		endfor
		
		############# use the weighted KMeans-like temporal component ###########################	
		case 1	
		############## UPDATE CENTROID TIMESTAMPS ########################
		#calculate the distance component for all observations
		if (deltaXmax ~= -1)
			mx = deltaXmax ^ 2;
		else
			mx = max(distObsClust);
		endif
		dists = 1 - gammaD * (distObsClust / mx);
	
		#first recompute current centroid timestamps
		for i=1:noClusters
			a =  (clusassign == i); 	#which individual belong in this cluster
			pointsInCluster(i) = sum(a);	#how many individuals are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
				fcts = dists(a);
				
				#calculate the obs for summing
				dt = reshape(dtimestamp(a), 1, length(dtimestamp(a)));
				selectedObs = dt .* fcts;	#numerator	
				selectedDiv = dt;		#put nans where there are nans in dt.
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* fcts;

				#calculate the new values			
				rez = sumskipnan( selectedObs );
				div = sumskipnan( selectedDiv );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				ctimestamp( i ) = rez / div;
			endif
		endfor

		############## UPDATE CENTROID POSITIONS ########################
		# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
		dts = reshape(dtimestamp, 1, length(dtimestamp)) .- oldCtimestamp(clusassign);
		fcts = dts .^ 2;

		if (deltaTmax ~= -1)
			mx = deltaTmax ^ 2;
		else
			mx = max(fcts);
		endif
		fcts = 1 - gammaT * (fcts ./ mx);
		
		for i=1:noClusters
			a =  (clusassign == i); 	#which individ belong in this cluster
			pointsInCluster(i) = sum(a);	#how many individ are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				#calculate the obs for summing
				selectedObs = data( a, :) .* repmat(fcts(a)', 1, columns(data));
				selectedDiv = data( a, :);
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* repmat(fcts(a)', 1, columns(data));

				#calculate the new values			
				rez = sumskipnan( selectedObs );
				div = sumskipnan( selectedDiv );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				centroids( i , : ) = rez ./ div;
			endif
		endfor
		
		############# When infering graph structure in clustering, the temporal update is dependent on the adjacency matrix ###########################
		case 3
		############## UPDATE CENTROID TIMESTAMPS ########################
		#calculate the multimensional distance component for all observations
		if (deltaXmax ~= -1)
			mx = deltaXmax ^ 2;
		else
			mx = max(distObsClust);
		endif
		dists = 1 - gammaD * (distObsClust / mx);
	
		for i=1:noClusters
			#which individual belong in this cluster, a vector as long as the observations, with one for those observations belonging to cluster i
			a =  (clusassign == i); 	
			pointsInCluster(i) = sum(a);	#how many individuals are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
				fcts = dists(a);
				
				# calculate the obs for summing. OBS: the following takes into account that there can be missing data in dtimestamp
				# missing data appears as NaN in the vector.
				dt = reshape(dtimestamp(a), 1, length(dtimestamp(a)));
				selectedObs = dt .* fcts;	# the numerator of the division
				selectedDiv = dt;		# put nans where there are nans in dt, so that they are not taken into account in the denominator
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* fcts;
				
				# calculate the temporal closeness component
				adjclus = Adj(:, i)';			# the adj(p, i) relation from all to current cluster
				adjclus = adjclus .^ 2;
				adjnum = oldCtimestamp .* adjclus; 	# part at the numerator
				adjclus(i) = [];	# remove components corresponding to current cluster i
				adjnum(i) = [];

				#calculate the new values			
				rez = lamb1 * sumskipnan( selectedObs ) + lamb2 * sumskipnan( adjnum );
				div = lamb1 * sumskipnan( selectedDiv ) + lamb2 * sumskipnan( adjclus );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				ctimestamp( i ) = rez / div;
			endif
		endfor	
		
		############## UPDATE CENTROID POSITIONS - identical typeT = 1 ########################
		# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
		dts = reshape(dtimestamp, 1, length(dtimestamp)) .- oldCtimestamp(clusassign);
		fcts = dts .^ 2;

		if (deltaTmax ~= -1)
			mx = deltaTmax ^ 2;
		else
			mx = max(fcts);
		endif
		fcts = 1 - gammaT * (fcts / (mx));
	
		for i=1:noClusters
			a =  (clusassign == i); 	#which individ belong in this cluster
			pointsInCluster(i) = sum(a);	#how many individ are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				#calculate the obs for summing
				selectedObs = data( a, :) .* repmat(fcts(a)', 1, columns(data));
				selectedDiv = data( a, :);
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* repmat(fcts(a)', 1, columns(data));

				#calculate the new values			
				rez = sumskipnan( selectedObs );
				div = sumskipnan( selectedDiv );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				centroids( i , : ) = rez ./ div;
			endif
		endfor
		
		###### For typeT = 4 algorithm, we use the temporal aware disimilarity distance between centroids, therefore complex centroid update ###############
		case {4,5}
		############## UPDATE CENTROID TIMESTAMPS ########################
		#calculate the multimensional distance component for all observations
		if (deltaXmax ~= -1)
			mx = deltaXmax ^ 2;
		else
			mx = max(distObsClust);
		endif
		dists = 1 - gammaD * (distObsClust / mx);
		
		# calculate the temporal distance from all centroids to all centroids
		ctdiff = distance(oldCentroids', oldCentroids') .^ 2;	# the SQUARRE of the Euclidean distance from each centroid to all the others
		
		if (deltaXmax ~= -1)
			divv = deltaXmax .^ 2;
		else
			divv = max(max(ctdiff));
		endif
		ctdiff = 1 - gammaD * (ctdiff ./ divv);
	
		# start recomputing timestamp for each centroid
		for i=1:noClusters
			#which individual belong in this cluster, a vector as long as the observations, with one for those observations belonging to cluster i
			a =  (clusassign == i); 	
			pointsInCluster(i) = sum(a);	#how many individuals are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
				fcts = dists(a);
				
				# calculate the obs for summing. OBS: the following takes into account that there can be missing data in dtimestamp
				# missing data appears as NaN in the vector.
				dt = reshape(dtimestamp(a), 1, length(dtimestamp(a)));
				selectedObs = dt .* fcts;	# the numerator of the division
				selectedDiv = dt;		# put nans where there are nans in dt, so that they are not taken into account in the denominator
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* fcts;
				
				# calculate the part corresponding to T2 (temporal-aware distance to other centroids)
				adjin = Adj(:, i)';			# the adj(p, i) relation from all to current cluster
				adjout = Adj(i, :);			# the adj(i, p) relation from current to all
				adjclus = ( adjin .^ 2 + adjout .^ 2 ) .* ctdiff(:, i)'; # this will also be at the denominator
				adjnum = oldCtimestamp .* adjclus ; 	# part at the numerator
				adjclus(i) = NaN;	# remove components corresponding to current cluster i
				adjnum(i) = NaN;

				#calculate the new values			
				rez = lamb1 * sumskipnan( selectedObs ) + lamb2 * sumskipnan( adjnum );
				div = lamb1 * sumskipnan( selectedDiv ) + lamb2 * sumskipnan( adjclus );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				ctimestamp( i ) = rez / div;
			endif
		endfor	
		
		############## UPDATE CENTROID POSITIONS ########################
		# calculate the temporal distance from all centroids to all centroids
		t1 = repmat( oldCtimestamp, length(oldCtimestamp), 1);
		t2 = repmat( oldCtimestamp', 1, length(oldCtimestamp));
		ctdiff = (t1 - t2) .^ 2;
		if (deltaTmax ~= -1)
			divv = deltaTmax .^ 2;
		else
			divv = max(max(ctdiff));
		endif
		ctdiff = 1 - gammaT * (ctdiff ./ divv);
		
		# calculate the temporal component (fcts(i) is the ponderation factor of obs i in the group)
		dts = reshape(dtimestamp, 1, length(dtimestamp)) .- oldCtimestamp(clusassign);
		fcts = dts .^ 2;

		if (deltaTmax ~= -1)
			mx = deltaTmax ^ 2;
		else
			mx = max(fcts);
		endif
		fcts = 1 - gammaT * (fcts / (mx));
	
		for i=1:noClusters
			a =  (clusassign == i); 	#which individ belong in this cluster
			pointsInCluster(i) = sum(a);	#how many individ are in this cluster
		
			if (pointsInCluster(i) > 0)	#if there are individuals assigned to this cluster, then update it
				#calculate the obs for summing
				selectedObs = data( a, :) .* repmat(fcts(a)', 1, columns(data));
				selectedDiv = data( a, :);
				selectedDiv( ~isnan(selectedDiv)) = 1;
				selectedDiv = selectedDiv .* repmat(fcts(a)', 1, columns(data));
				
				#calculate the term corresponding to T2 (the temporal-aware distance to other centroids)
				adjin = Adj(:, i)';			# the adj(p, i) relation from all to current cluster
				adjout = Adj(i, :);			# the adj(i, p) relation from current to all
				adjclus = ( adjin .^ 2 + adjout .^ 2 ) .* ctdiff(:, i)'; # this will also be at the denominator
				adjclus(i) = [];				#remove component correponding to the current cluster
				# right now, adjclus is a vector corresponding to each cluster except the current, therefore dimension noClusters-1
				# each clusters contributes to the current cluster proportional to its correponding weight in adjclus
				# need to replicate adjclus on each dimension in order to weight each dimension
				weights = repmat(adjclus', 1, columns(data));	# matrix of dimension (noClusters-1, columns(data))
				
				# select all clusters, except current one
				adjnum = oldCentroids;
				adjnum(i, :) = [];
				
				adjnum = adjnum .* weights ; 	# part at the numerator

				#calculate the new values			
				rez = lamb1 * sumskipnan( selectedObs ) + lamb2 * sumskipnan( adjnum ) ; # a row vector with length = columns(data), no features
				div = lamb1 * sumskipnan( selectedDiv ) + lamb2 * sumskipnan( weights );
				div(div == 0) = eps;	#to avoid 0/0 when all members of a column are NaNs
				centroids( i , : ) = rez ./ div;
			endif
		endfor
	
		############# ALL THE OTHER CASES ###########################	
		otherwise
			error(sprintf("I don't recognize 'typeT=%d' Did you give me a wrong typeT algorithm to work with???", typeT));
	endswitch
endfunction
