## Based on the drawGraph.m file fromt the PMTK3 toolbox
## drawGraph Automatic graph layout: interface to Neato (see http://www.graphviz.org/)
##
## Mandatory arguments:
##	inputfile - files where data usefull for the matrix construction will be read from (will be read with "loadCsvDataMissing")
##			OR
##	files - the files part when read with "loadCsvDataMissing"
##	tags - the tags part
##	vals - the values part
##
##	outputfile - the name of the DOT output file
##
##
## Optional arguments (string/value pair) [default in brackets]
##	'graph_threshold' [ -1 ] - filters the constructed Adjacency matrix for values inferior to threshold. Default -1 (disabled), in which 
##				case, the value of Adj(1,1) is used.
##	'omega1' and 'omega2' - are the weights of the temporal distance and the intersection distance between clusters in
##				the Adjacency matrix. Passed to the "createAdjacencyMatrix" function.
##	'omega' [ -1 , off] - a single parameter to set the values of omega1 = 1 - omega and omega2 = omega
## 	'labels' - labels{i} is a *string* for node i [1:n]
## 	'removeSelfLoops' - [1]
## 	'removeIsolatedNodes' - [1]
## 	'directed' - [default: determine from symmetry of Adj] 
## 	'systemFolder' - [defaults to location returned by graphvizRoot()]
##	'useDayMonthYear' [ false ] - transform the cluster timestamp in dd/mm/yy format?
##	'ctimestamp' [ [], empty ] - the centroids timestamp. If empty, they will be recomputed as the averages of timestamps of observations in clusters
##	'centroidsFile" [ [] ] - gives the file in which the centroids are stored. Used to load the centroid timestamp. Usefull in conjunction with 'inputfile'
##	'outputGephi' [false] - if true, output the 'gephi.nodes.csv' and 'gephi.edges.csv' to be imported into Gephi (www.gephi.org).
##
## Example
##	generateGraphFile('inputfile', 'output-clustered.csv', 'outputfile', 'output-graph.dot');
##	[f,t,v] = loadCsvFileMissing('output-clustered.csv'); generateGraphFile('files', f, 'tags', t, 'vals', v, 'outputfile', 'output-graph.dot');	
##
##	generateGraphFile('inputfile', 'output-clustered.csv', 'outputfile', 'output-graph.dot', 'labels', {'a','bbbb','c','d','e'})
##	generateGraphFile('inputfile', 'output-clustered.csv', 'outputfile', 'output-graph.dot',  'removeIsolatedNodes', 0, 'removeSelfLoops', 0);
##
## Modified:	29/09/2013 Improved the help, added receiving and passing the ctimestamp
##		13/11/2013 Added exportation to Gephi
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [graph_threshold, no_arcs] = generateGraphFile(varargin)

if ~exist('graphvizRoot','file')
  systemFolder = [];
  str = sprintf('%s %s %s\n', ...
		'warning: you should put the file graphvizRoot.m', ...
		'into the graphviz/bin directory', ...
		'if you  want to convert to postscript');
  %fprintf(str)
else
  %'C:\Program Files\ATT\Graphviz\bin', ...
  systemFolder = graphvizRoot();
end

[systemFolder, labels, removeSelfLoops, Adj, removeIsolatedNodes, directed, inputfile, outputfile, files, tags, vals, edgedescription, useDiameters, useDayMonthYear, filtcycle, graph_threshold, omega1, omega2, omega, ctimestamp, outputGephi, centroidsFile, no_arcs] = ...
    process_options(varargin, 'systemFolder', systemFolder, ...
		'labels', [], 'removeSelfLoops', 1, 'Adj', [], 'removeIsolatedNodes', false, 'directed', [], ...
		'inputfile', [], 'outputfile', [], 'files', [], 'tags', [], 'vals', [], 'edgedescription', true , ...
		'useDiameters', true, 'useDayMonthYear', false, 'filtcycle', false, 'graph_threshold', -1, 'omega1', 1, ...
		'omega2', 1, 'omega', -1, 'ctimestamp', [], 'outputGephi', false, 'centroidsFile', [], 'no_arcs', -1 );

if (isempty(inputfile) && isempty(files) && isempty(tags) && isempty(vals) )
	
	error('You need to give me at least the file from where to read my data!');
endif

if (isempty(outputfile))
	outputfile = 'output-graph.dot';
	warning(sprintf("No DOT output file given! Outputing to file \"%s\"", outputfile));
endif

# check if omega was set
if (omega > -1)
	omega1 = 1 - omega;
	omega2 = omega;
endif

if (outputGephi)
	gephiNodes = fopen("gephi.nodes.csv", "w");
	fprintf(gephiNodes, "id,timestamp\n");
	gephiEdges = fopen("gephi.edges.csv", "w");
	fprintf(gephiEdges, "Source,Target,noEntities,weight\n");
endif

########################## DONE CONFIGURING SYSTEM STUFF ########################

# first we need to construct the Adjacency matrix
if ( ~isempty(inputfile) )
	printf("--> Loading data from file '%s' ... ", inputfile);
	[files,tags,vals] = loadCsvFileMissing( inputfile, true );
	printf(" done!\n");
endif

# last column is the cluster assignement, files is the file assignement, dates are the first column
dtimestamp = vals(:, 1);
clusassign = vals(:, end);

# from the centroids file we read only their timestamp
if (~isempty(centroidsFile))
	if ( ~isempty(ctimestamp) )
		warning("Both 'ctimestamp' and 'centroidsFile' were set. Ignoring 'centroidsFile'");
	else
		# load centroid timetamp from file
		[foo1, foo2, cval] = loadCsvFileMissing( centroidsFile );
		ctimestamp = cval(:, 1);
	endif
endif

# If the cluster timestamp was not given, it is recalculated here.
if ( isempty(ctimestamp) )
	ctimestamp = zeros(1, max(clusassign));
	for i = 1:max(clusassign)
		pos = find(clusassign == i);
		ctimestamp(i) = sum(dtimestamp(pos)) / length(pos);
	endfor
endif

# if it is not given, calculate the Adjacency matrix, using the post clustering method proposed in IJAIT 14.
if (isempty(Adj))
	[Adj, Adj_count, arcs, countr, startsat] = createAdjacencyMatrix('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', filtcycle, 'omega1', omega1, 'omega2', omega2, 'useDiameters', useDiameters, 'ctimestamp', ctimestamp);
else
	[foo1, Adj_count, foo2, arcs, entities, startsat] = interPhi('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', filtcycle, 'ctimestamp', ctimestamp);
#	Adj_count
#	interphi = foo1
endif

# calculate the threshold to filter the Adjacency matrix
if ((graph_threshold == -1 ) && (no_arcs == -1) )#none given
	# one option is to get the Adj(1,1) as the threshold. This ensures that there are no self-arcs
#	graph_threshold = Adj(1,1);

	# or we can use the heuristic to choose the highest n-1 arcs, where n is the number of nodes
	no_arcs = rows(Adj) - 1;
endif

# if the threshold was not specified, then we select no_arcs arcs from the matrix, those with highest scores
# to do this, we compute the threshold needed to filter out the other arcs
if ((graph_threshold == -1 ) && (no_arcs != -1) )
	Adj_vals = sort(reshape(Adj, 1, []), 'descend') ; # all the values in the adjacency matrix, sorted inversly
	
	# the value of the threshold becomes that of the (no_arcs+1)th arc, since the filtering is inclusive
	graph_threshold = Adj_vals(no_arcs + 1);
	printf("Threshold automatically calculated for %d arcs: %.5f\n", no_arcs, graph_threshold);
	
	Adj_vals(Adj_vals <= graph_threshold) = 0;
	if (~isempty(find(Adj_vals)))
		no_arcs = find(Adj_vals)(end);
	endif
endif


################### transform the Adjacency matrix ###############################

# remove self loops
if removeSelfLoops
  Adj = setdiag(Adj, 0); 
end

[n,m] = size(Adj);
if n ~= m, warning('not a square Adjacency matrix!'); end
%if ~isequal(diag(Adj),zeros(n,1)), warning('Self-loops in Adjacency matrix!');end

if isempty(labels)
  labels = cell(1,n);
  for i=1:n
    labels{i} = sprintf('%d', i); 
  end
end

# determine if it is a directed graph
if isempty(directed)
  %if isequal(triu(Adj ,1),tril(Adj,-1)'), directed = 0; else, directed = 1; end 
  if isequal(Adj,Adj')
    directed = 0;
  else
    directed = 1;
  end
end

# filter the calculated Adjacency matrix
Adj_score = Adj .* 100;	# save the initial Adjacency matrix TODO multiplied the score by 100 to see it on the matrix
Adj( Adj <= graph_threshold ) = 0;	# filter the stuff scores under threshold

# remove isolated nodes
if removeIsolatedNodes
  isolated = [];
  for i=1:n
    nbrs = [find(Adj(i,:)) find(Adj(:,i))'];
    nbrs = setdiff(nbrs, i);
    if isempty(nbrs)
      isolated = [isolated i];
    end
  end
  Adj = removeRowsCols(Adj, isolated, isolated);
  Adj_count = removeRowsCols(Adj_count, isolated, isolated);
  arcs = removeRowsCols(arcs, isolated, isolated);
  labels(isolated) = [];
end

# I tried to filter arcs so that only those having entities remain.
# but it boils down to setting omega1 = 0
#select = Adj_count > 0;
#Adj = Adj .* select;

Adj = double(Adj > 0); 	# make sure it is a binary matrix: cast to double type

########## end of transformation of Adjacency matrix ###############

# calculate ranks
ranks = rankGraph(Adj);

#here started the graph_to_dot function
#set the new parameters from the old function
node_label = labels;   arc_label = arcs;
width = 41;        height = 31;
leftright = 1;
showcountries = 0; #ctimestamp = [];
########## end of configuration ###############

fid = fopen(outputfile, 'w');
if fid==-1
  error(sprintf('could not write to %s', outputfile))
end
if directed
    fprintf(fid, 'digraph G {\n');
    arctxt = '->'; 
    if isempty(arc_label)
        labeltxt = '';
    else
        labeltxt = '[label="%s", weight=1 ]';
    end
else
    fprintf(fid, 'graph G {\n');
    arctxt = '--'; 
    labeltxt = '[dir=none]';
    if ( ~isempty(arc_label))
        labeltext = '[label="%s",dir=none, weight=1 ]';
    end
end

fprintf(fid, 'center = 1;\n');
fprintf(fid, 'size=\"%d,%d\";\n', width, height);
if leftright
    fprintf(fid, 'rankdir=LR;\n');
end
#fprintf(fid, 'ranksep=4.0;\n');
fprintf(fid, 'nodesep=1.0;\n');

Nnds = length(Adj);

if (showcountries)
	#put the countries in the graph
	fprintf(fid, '{\n\tnode [shape=box];\n\t{rank = same; ');
	for i=1:length(countr)
		fprintf(fid, '"%s"; ', countr{i});
	endfor
	fprintf(fid, '}\n}\n\n');

	#and the links of the countries
	for i=1:length(countr)
		fprintf(fid, '"%s" -> %d [color="cyan1"]; \n', countr{i}, startsat(i));
	endfor
	fprintf(fid, '\n');
endif

Nnds = length(Adj);

# the color for the nodes
	colors={"blue", "cyan", "red", "magenta", "green", "yellow", "cyan", "blue", "black", "red", "magenta", "green", "yellow", "black", "blue", "cyan", "red", "magenta", "green", "yellow", "black", "blue", "cyan", "red", "magenta", "green", "yellow", "black"}; #switched the second cyan with black
	markers={"+", "*", "o", "x", "^", "+", "*", "o", "*", "o", "+", "^", "*", "o","+", "*", "o", "x", "^", "+", "*", "o", "*", "o", "+", "^", "*", "o"};
for node = 1:Nnds               % process NODEs 
    extra = "";
    color = sprintf("color=\"%s\", style=filled", colors{node});
    if (node == 9) # switched this for the european companies dataset
        color = [ color sprintf(", fontcolor = white")];
    endif
    
    if isempty(node_label)
        fprintf(fid, '%d [ %s ] ; %s\n', node, color, extra);
    else
        fprintf(fid, '%d [ label = <<b>%s | %s</b>>, %s ]; %s\n', node, node_label{node}, printDate(ctimestamp(node), useDayMonthYear), color, extra);
    endif
    
    if (outputGephi)
    	fprintf(gephiNodes, "%d,%s\n", node, printDate(ctimestamp(node), useDayMonthYear));
    endif
end

edgeformat = strcat(['%d ',arctxt,' %d ',labeltxt,';\n']);
for node1 = 1:Nnds              % process ARCs
    if directed
        arcs = find(Adj(node1,:));         % children(Adj, node);
    else
        arcs = find(Adj(node1,node1+1:Nnds)) + node1; % remove duplicate arcs
    end
    for node2 = arcs
	if isempty(arc_label)
		fprintf(fid, edgeformat, node1, node2);
	else
		if (edgedescription)
			# if "edgedescription" is set then print on the edge the list of entitites
	        	fprintf(fid, edgeformat, node1, node2, arc_label{node1, node2} );
	        else
	        	# if not, print how many individuals are passing from state i to state j
	        	fprintf(fid, edgeformat, node1, node2, sprintf("%d entities | score %.2f", Adj_count(node1, node2), Adj_score(node1, node2)) );	        	
	        endif
        endif
        
    	if (outputGephi)
    		fprintf(gephiEdges, "%d,%d,%d,%.4f\n", node1, node2, Adj_count(node1, node2), Adj_score(node1, node2) );
   	 endif
    endfor
end
fprintf(fid, '}'); 
fclose(fid);

################## done generating the dot file, now transform to PNG and PDF

[DIR, NAME, foo, foo] = fileparts (outputfile);
if (length(DIR) == 0)	# handle the case when file is in the current dir
	DIR = ".";
endif
NAME = [DIR '/' NAME]; #include in the name also the path. Now from "/path/output.dot", NAME will be "/path/output"

#transform to PNG
cmd = sprintf('dot -Tpng -o "%s.png" "%s.dot"', NAME, NAME);
status = system(cmd);
if status ~= 0
	error(sprintf('error executing %s', cmd));
end
  
#transform to PDF, it was initially PS
dot_to_pdf(outputfile, sprintf('%s.pdf', NAME), false); #pdf doesn't like the landscape option
  
#transform te PS to PDF
if ~isempty(outputfile)
#	if ismac
#		str = '/usr/local/bin/ps2pdf';
#	else
#		str = 'ps2pdf';
#	endif
#	cmd = sprintf('%s %s.ps %s.pdf', str, NAME, NAME);
#	status = system(cmd);
#	if status ~= 0
#		error(sprintf('error executing %s', cmd));
#	endif
#	  
#	  
#	cmd = sprintf('rm %s.ps', NAME);
#	status = system(cmd);
#	if status ~= 0
#		error(sprintf('error executing %s', cmd));
#	endif

	cmd = sprintf('pdfcrop "%s.pdf" > /dev/null && mv "%s-crop.pdf" "%s.pdf"', NAME, NAME, NAME);
	status = system(cmd);
	if status ~= 0
		error(sprintf('error executing %s', cmd));
	endif
 
end

# close gephi files if demanded.
if (outputGephi)
	fclose(gephiNodes);
	fclose(gephiEdges);
endif
endfunction #generateGraphFiles

%%%%%%%%%%

function M = removeRowsCols(M, rows, cols)
% Remove rows and columns from a matrix
% Example
% M = reshape(1:25,[5 5])
%> removeRowsCols(M, [2 3], 4)
%ans =
%     1     6    11    21
%     4     9    14    24
%     5    10    15    25
     
[nr nc] = size(M);

ndx = [];
for i=1:length(rows)
  tmp = repmat(rows(i), nc, 1);
  tmp2 = [tmp (1:nc)'];
  ndx = [ndx; tmp2];
end
for i=1:length(cols)
  tmp = repmat(cols(i), nr, 1);
  tmp2 = [(1:nr)' tmp];
  ndx = [ndx; tmp2];
end
if isempty(ndx), return; end
k = subv2ind([nr nc], ndx);
M(k) = [];
M = reshape(M, [nr-length(rows) nc-length(cols)]);

%%%%%%%%%
end

function M = setdiag(M, v)
% SETDIAG Set the diagonal of a matrix to a specified scalar/vector.
% M = set_diag(M, v)

n = length(M);
if length(v)==1
  v = repmat(v, 1, n);
end

% e.g., for 3x3 matrix,  elements are numbered
% 1 4 7 
% 2 5 8 
% 3 6 9
% so diagnoal = [1 5 9]


J = 1:n+1:n^2;
M(J) = v;

%M = triu(M,1) + tril(M,-1) + diag(v);
end

%%%%%%%%%%%

function ndx = subv2ind(siz, subv)
% SUBV2IND Like the built-in sub2ind, but the subscripts are given as row vectors.
% ind = subv2ind(siz,subv)
%
% siz can be a row or column vector of size d.
% subv should be a collection of N row vectors of size d.
% ind will be of size N * 1.
%
% Example:
% subv = [1 1 1;
%         2 1 1;
%         ...
%         2 2 2];
% subv2ind([2 2 2], subv) returns [1 2 ... 8]'
% i.e., the leftmost digit toggles fastest.
%
% See also IND2SUBV.

 
if isempty(subv)
  ndx = [];
  return;
end

if isempty(siz)
  ndx = 1;
  return;
end

[ncases ndims] = size(subv);

%if length(siz) ~= ndims
%  error('length of subscript vector and sizes must be equal');
%end

if all(siz==2)
  %rbits = subv(:,end:-1:1)-1; % read from right to left, convert to 0s/1s
  %ndx = bitv2dec(rbits)+1; 
  twos = pow2(0:ndims-1);
  ndx = ((subv-1) * twos(:)) + 1;
  %ndx = sum((subv-1) .* twos(ones(ncases,1), :), 2) + 1; % equivalent to matrix * vector
  %ndx = sum((subv-1) .* repmat(twos, ncases, 1), 2) + 1; % much slower than ones
  %ndx = ndx(:)';
else
  %siz = siz(:)';
  cp = [1 cumprod(siz(1:end-1))]';
  %ndx = ones(ncases, 1);
  %for i = 1:ndims
  %  ndx = ndx + (subv(:,i)-1)*cp(i);
  %end
  ndx = (subv-1)*cp + 1;
end
end
%%%%%%%%%%%

function d = bitv2dec(bits)
% BITV2DEC Convert a bit vector to a decimal integer
% d = butv2dec(bits)
%
% This is just like the built-in bin2dec, except the argument is a vector, not a string.
% If bits is an array, each row will be converted.

[m n] = size(bits);
twos = pow2(n-1:-1:0);
d = sum(bits .* twos(ones(m,1),:),2);

end

function status = dot_to_pdf(dotname, outname, landscape)

               
% Useful options:
%   -Glandscape (outputs in landscape mode)
%   -Gconcentrate (merges two-way edges into one way edge, displays
%   parallel edges in different way)
%   -Gratio=.707 (changes to A4 landscape aspect ratio, other options are "fill", "compress",
%   "expand", "auto")
%   -Ncolor="blue" (changes node outlines to blue)
%   -Ecolor="red" (changes edges to red)
%   -Earrowsize=2 (changes size of arrows)
%   -Nstyle="filled" -Nfillcolor="#ddddff" (make nodes light blue)
%   -Nfontsize=32 (change font size to 32pt)
%   -Gnodesep=0.125 (make nodes twice as close
%   -Nshape="box" (change node shape to box)
%
% Details here:
% http://www.graphviz.org/doc/info/attrs.html

%opts = ' -Gconcentrate -Gratio=.707 -Ncolor="blue" -Ecolor="green" -Earrowsize=2 -Nstyle="filled" -Nfillcolor="#ddddff" -Nfontsize=40 ';
#opts = '-Gconcentrate  -Gratio=0.707 -Earrowsize=2 -Nfontsize=25 -Epenwidth=3';
opts = '-Gconcentrate';

if landscape
    opts = strcat(opts,' -Glandscape ');
end

if ismac
  dotstr = '/usr/local/bin/dot';
else
  dotstr = 'dot';
end
%cmd = strcat('C:\temp\graphviz-2.8\bin\dot ',opts,' -T ps -o graphVizIt.ps graphVizIt.txt ')
cmd = sprintf('%s %s -Tpdf "%s" -o "%s"', dotstr, opts, dotname, outname);
%cmd = sprintf('dot -Tps %s -o %s', dotname, outname);
%cmd = sprintf('dot %s -Tpng %s -o %s', opts, dotname, outname);
status = system(cmd);
if status ~= 0
  error(sprintf('error executing %s', cmd));
end

end


# function for printing date, 
function datestring = printDate(ctimestmp, useDayMonthYear)
	if (useDayMonthYear)
		datestring = datestr(ctimestmp, "dd/mm/yy");
	else
		datestring = sprintf("%d", ctimestmp);
	endif
endfunction
