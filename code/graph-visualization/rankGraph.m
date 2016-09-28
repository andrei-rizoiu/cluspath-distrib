function [ranks] = rankGraph(adj)

	dim = rows(adj);
	ranks = zeros(1, dim);
	father = zeros(1, dim);
	queue = [];
	
	indeg = sum( adj > 0 );
	queue = find( indeg == min( indeg ) );
	ranks(queue) = 1;
	
	while ( ~isempty(queue) )
		elem = queue(1);
		queue(1) = [];
		
		for suc=1:dim
			if (adj(elem, suc) ~= 0)	#if suc is a succesor of elem
			
				if (ranks(suc) == 0)	#not seen before
					ranks(suc) = ranks(elem) + 1;
					queue = [queue suc];
					father(suc) = elem;
				endif
			
				if ( (ranks(suc) ~=0 ) && (ranks(suc) <= ranks(elem)) )
					#check out if suc is a predecesor of elem
					fath = father(elem);
					while ( (fath ~= suc) && (fath ~= 0) )
						fath = father(fath);
					endwhile
					
					if fath == 0	#suc is not a predecesor of elem
						ranks(suc) = ranks(elem) + 1;
						queue = [queue suc];
						father(suc) = elem;
					else		#suc is a predecesor of elem
						#put all on the same cycle at the same rank
						fath = father(elem);
						while ( fath ~= suc )
							ranks(fath) = ranks(suc);
							fath = father(fath);
						endwhile
						ranks(elem) = ranks(suc);
					endif
				endif
			endif
		endfor
	endwhile
endfunction
