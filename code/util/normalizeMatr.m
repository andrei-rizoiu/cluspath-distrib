## Usage: [uN] = normalizeMatr(V, standardize)
##
## Return the the normalized vectors in a matrix. All values are put in [0, 1].
## Each column is considered to be a vector.
## If standardize is true (default) then all the columns of V are standardized
## x(i) = ( x(i) - aver(x) ) / standard_deviation(x)

function [uN] = normalizeMatr(V, standardize)

	if (nargin < 2)
		standardize = true;
	endif

	if ( standardize )
		meanNH = repmat(mean(V), rows(V), 1);
		stdN = std(V);
		stdN(stdN == 0) = eps;
		stdNH = repmat(stdN, rows(V), 1);
		uN = (V - meanNH ) ./ stdNH;
	else
		maxN = max(V);
		minN = min(V);
		widthN = maxN - minN;
		#there is a possibility that the width is 0 (when all the elements are equal or NaN)
		#replace the 0 by eps and let the division make the choice
		widthN(widthN == 0) = eps;
	
		minNH = repmat(minN, rows(V), 1);
		widthNH = repmat(widthN, rows(V), 1);
	
		uN = (V - minNH) ./ widthNH;
	endif
endfunction
