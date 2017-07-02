%load edge list of newqtork - textfile in the format 
EdgeList = load('facebook_combined.txt'); 

%sort the edge list
sortrows(EdgeList);

% number of edges in the graph
NUMBER_EDGES = length(EdgeList);

% get all nodes, sorted
nodes = sort(unique([EdgeList(:, 1) EdgeList(:, 2)]));

% number of nodes in the graph
NUMBER_NODES = numel(nodes);

% initialize adjacency matrix
Adj = zeros(NUMBER_NODES);
    
% iterate across all edges and populate adjacency matrix with one's
for d = 1 : NUMBER_EDGES; 
    Adj(find(nodes == EdgeList(d, 1)), find(nodes == EdgeList(d, 2))) = 1; 
    Adj(find(nodes == EdgeList(d, 2)), find(nodes == EdgeList(d, 1))) = 1;
end    

%make adjacency matrix sparse
Adj = sparse(Adj);

% Finally, the adjacency matrix should be 1 whereever an edge exits, but
% the construction will have greated entries greater than one.
% This line of code, firstly takes all the non-zeros in A and sets them
% to logical one: (A>0) does this.  Finally the multiplication by 1.0 makes
% the non-zeros into floating-point ones, which is required by the eigs 
% function later on.
Adj = 1.0*(Adj>0);


%Form the (normalised) Laplacian matrix
% d is the vector of vertex degrees
d = sum(Adj);

% D is a diagonal matrix whose diagonal is set to the degrees (note, we
% keep D sparse).
D = sparse(1:NUMBER_NODES, 1:NUMBER_NODES, d);

% Dinv in the inverse of D
Dinv = sparse(1:NUMBER_NODES,1:NUMBER_NODES, 1./d);

% L is the Laplacian matrix
L = D - Adj;

% nL is the normalised Laplacian matrix
nL = Dinv*(D - Adj);

% Use the eigensolver for sparse matrices, to calculate the
% 3 eigenvectors of smallest magnitude
[V, E] = eigs(nL, 3, 'sm');

% Use the second and third eigenvectors to form co-ordinates
% and plot the resulting graph using gplot.
W = V(:,2:3);

%index the vector starting from 1 instead 0
%for index of the vector here is used the id's of the vertices,
%which in this case start from 0
for d = 1 : NUMBER_NODES; 
  nodes(d) = nodes(d) + 1;
end

NUMBER_K_COMMUNITIES = 9;

fid = fopen('results.txt', 'w');
header1 = 'Community number ';
header2 = 'Modularity';
header3 = 'Number of edges cut';
fprintf(fid, [header1  '   ' header2 '   ' header3 '\n']);
%compute 2 - k communities
communities = {};
  
for k = 2 : NUMBER_K_COMMUNITIES; 
  % Use kmeans to cluster the points in the vector W
  %[classes, centers, sumd, D] = kmeans(W,4);
  clust = kmeans(W, k);
  

  for kk = 1 : k; 
    %for each community - create a vector with vertices in it
    communities{kk} = nodes(clust==kk, 1);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Compute modularity per each  k kmeans cutting %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Q = 0;
  communityEdges = 0;
  for m = 1 : length(communities);
    communityEdges = sum(sum(Adj(communities{m}, communities{m})))/2; 
    e_mm = communityEdges / NUMBER_EDGES;
    a_m = sum(sum(Adj(:, communities{m}))) / NUMBER_EDGES - e_mm;
    Q = Q + (e_mm - a_m^2); 
  end
  
  % number of edge cuts per community NUMBER_EDGES-communityEdges
  fprintf(fid, '%f             %f              %f \n', [k Q NUMBER_EDGES-communityEdges]');
end

fclose(fid);









