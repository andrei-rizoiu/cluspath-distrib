function [d] = calculateDiameter(data, missing)

	if (nargin < 2)
		missing = false;
	endif

	if ( ! missing )
		dist = distance(data', data');
		if ( sum(sum(isnan(data))) > 0)
			#there are NaNs in the distance, maybe someone forgot to tell us that we have missing data

#			warning("You forgot to tell us that there is missing data in the dataset!");
			missing = true;
			dist = distanceMissing(data', data');
		endif
	else
		dist = distanceMissing(data', data');
	endif

	d = max(max(dist));
endfunction
