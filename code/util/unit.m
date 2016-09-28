## Usage: [uV] = unit(V)
##
## Return the unit vectors in a matrix. They keep the direction and modulo 1.
## Each column is considered to be a vector

function [uV] = unit(V)
	[foo1, foo2, SSQ] = sumskipnan(V);
	normV = sqrt(SSQ);
	#there is a possibility that the norm is 0 (when all the elements are either 0 or NaN)
	#replace the 0 by eps and let the division make the choice
	normV( normV == 0) = eps;
	normVH = repmat(normV, rows(V), 1);
	uV = V ./ normVH;
endfunction
