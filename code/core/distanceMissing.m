## Usage: d = distanceMissing(a,b)
##
## DISTANCE - computes Euclidean distance matrix when there are missing values (NaN).
## OBS: there will be serious speed penalisations compared to the distance(a, b) function.

function d = distanceMissing(a,b)

if (nargin ~= 2)
   error('Not enough input arguments');
end

if (size(a,1) ~= size(b,1))
   error('A and B should be of same dimensionality');
end

d = [];
for i=1:size(a, 2)
	da = repmat( a(:, i), 1, size(b,2));
	
	dd = sqrt( sumskipnan((da-b) .* (da-b),1));
	d = [d ; dd];
endfor

