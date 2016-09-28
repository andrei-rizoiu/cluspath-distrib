## Just a wrapper for the temporalKMeans function to be used with the parcellfun in genetic_pareto_front

function [Adj, clusassign, centroids, ctimestamp, centroids_EVOLUTION] = temporalKMeans_wrapper(params, centroids)

#	Adj = params.Alpha;
#	clusassign = params.Beta;
#	centroids = params.Delta;
#	ctimestamp = zeros(10);
	[Adj, clusassign, dtimestamp, centroids, ctimestamp, centroids_EVOLUTION] = temporalKMeans( ...
        			'file', params.dataset, 'typeT', params.typeT, 'noClusters', params.noClusters, 'maxiter', params.maxiter, ...
        			'Alpha', params.Alpha, 'Beta', params.Beta, 'Delta', params.Delta, ...
        			'lamb1', params.lamb1, 'lamb2', params.lamb2, 'lamb3', params.lamb3, 'loadCentroidsFromFile', centroids, ...
        			'saveResults', false, 'outputfid', 0, 'saveInitCentroids', false); #, ... 'useDiameters', false 
	
endfunction
