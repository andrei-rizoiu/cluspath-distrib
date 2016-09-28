## Just a wrapper for the evaluateTemporalClustering function to be used with the parcellfun in genetic_pareto_front

function [mdvar, tempvar, shannonP, smooth] = evaluateTemporalClustering_wrapper(params, Adj, clusassign, centroids, ctimestamp)

	[mdvar, tempvar, shannonP, smooth] = evaluateTemporalClustering('inputfile', params.dataset, 'Alpha', params.Alpha, ...
		'clusassign', clusassign, 'centroids', centroids, 'ctimestamp', ctimestamp, 'Adj', Adj, ...
		'outputfid', 0 );

endfunction
