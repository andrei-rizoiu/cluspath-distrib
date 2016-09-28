## Usage: [uN] = deNormalizeMatr(u, origV, standardize)
##
## Returns the the denormalized vectors in a matrix. Values are put from [0, 1] to an unnormalised form.
## Each column is considered to be a vector. The de-normalization is done in the space defined by origV. 
## This is used to return the dates of centroids from the normalised space to the original space.
## ex: origX = X; XN = normalizeMatr(X); XDN = deNormalizeMatr(XN, X);
## XDN will be the matrix denormalized in the space defined by X. Then XDN must be equal to X
##
## If standardize is true (default) then all the columns of V are standardized
## x(i) = ( x(i) - aver(x) ) / standard_deviation(x)
## This previous formula is undone.

function [uN] = deNormalizeMatr(u, origV, standardize)

	# check that u and origV have the same number of dimensions
	if(size(u, 2) ~= size(origV, 2) )
		error("u is not defined in the same space as origV! number of columns differ %d ~= %d!", size(u, 2), size(origV, 2));
	endif

	# default to standardization
	if (nargin < 3)
		standardize = true;
	endif

	if ( standardize )
		meanNH = repmat(mean(origV), rows(u), 1);
		stdN = std(origV);
#		stdN(stdN == 0) = eps; # not applicable for destandardization
		stdNH = repmat(stdN, rows(u), 1);
		uN = u .* stdNH + meanNH ;
	else
		maxN = max(origV);
		minN = min(origV);
		widthN = maxN - minN;
		#there is a possibility that the width is 0 (when all the elements are equal or NaN)
		#replace the 0 by eps and let the division make the choice
#		widthN(widthN == 0) = eps; # not applicable for denormalization
	
		minNH = repmat(minN, rows(u), 1);
		widthNH = repmat(widthN, rows(u), 1);
	
		uN = u .* widthNH + minNH ;
	endif
endfunction
