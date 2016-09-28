## Function to calculate if the domination count (how many other individuals dominate a given individual)
## for a set of individuals. All criteria (given by the number of rows in the matrix X) need to be minimised.
## It is not meant to becalled directly, but rather from the "genetic_pareto_front" function

function f = domination_count(X, optimizationCriteria)

f = zeros(1, columns(X));
for i=1:columns(X)
    for j=1:columns(X)
        if ( (i ~= j) && (prod(X(optimizationCriteria, i) <= X(optimizationCriteria, j)) * max(X(optimizationCriteria, i) < X(optimizationCriteria, j)) == 1) )
		f(j) = f(j) + 1;
        end
    end
end

