function [ newcentroids, newpos] = eliminateDuplicates( centroids )

# the function receives a set of initial centroids and it verifies that they are all distinct
# if not prunes the centroids and returns the distinct list and their positions	

newcentroids = centroids ;
newpos = 1:rows(centroids);

E = distance( centroids', centroids' );
# fill up the diagonal
E = E + tril(ones( rows(E) ) ) ;
zerouri = sum(E == 0 );
if ( sum(zerouri > 0) ~=0 )
	# means we have at least duplicate
	printf("Duplicate initial centroid detected... correcting\n");
	pos = (zerouri == 0);
	newcentroids = centroids ( pos, :);
	newpos = 1:rows(centroids);
	newpos = newpos(pos);
endif

end
