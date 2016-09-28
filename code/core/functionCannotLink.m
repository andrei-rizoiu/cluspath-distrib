## Usage: [CLP] = functionCannotLink(tdata, clusassign, args)
##
## Function CannotLinkPenality calculats the penality based on the
## cannotlink component. CLP(i,j) is the penality if the element di is in cluster cj
## The function is: alpha * |t(di) - t(dj)|
## Complete function: alpha * sum(dk) [|t(di) - t(dk)| * (clus(di) == clus(dk)) ]

function [CLP] = functionCannotLink(tdata, clusassign, noclusters, args)
	
	if (nargin < 3)
		args = [];
	endif
	
	if ( isempty(args) )
		alpha = 1;
	else
		alpha = args(1);
	endif
		
	noindivid = length(tdata);
	CLP = zeros(noindivid, noclusters);
	if (alpha == 0)
		return;
	endif
	
	## to calculate CLP(i,j), we must find the individuals dk so that clus(di) == clus(dk)
	## then calculate sum(dk)|t(di) - t(dk)|
	## multiply by alpha
	
	for i=1:noindivid
		temp = repmat(vec(clusassign), 1, noclusters);
		temp(i, :) = [1:noclusters];
		temp = (temp == repmat(temp(i, :), noindivid, 1)) + 1 - 1;
		temp( temp == 0) = NaN;
		temp = temp .* repmat(vec(tdata), 1, noclusters);
		temp = abs(tdata(i) - temp);	#here the part of individual i gets to be zero
		[temp, Ni] = sumskipnan(temp);
		if ( (Ni-1) > 0 )
			#one of the lines is exactly the individual i, which contributes with zero, but is counted in Ni
			temp = temp ./ (Ni - 1);
		endif
		CLP(i, :) = temp;
	endfor
	CLP = CLP * alpha;
	
endfunction
