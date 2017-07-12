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

% Use the second and third eigenvectors 
W = V(:,2:3);

%index the vector starting from 1 instead 0
%for index of the vector here is used the id's of the vertices,
%which in this case start from 0
for d = 1 : NUMBER_NODES; 
  nodes(d) = nodes(d) + 1;
end

NUMBER_K_COMMUNITIES = 3;

fid = fopen('result.txt', 'w');
header1 = 'Community number ';
header2 = 'Modularity';
header3 = 'Number of edges cuts';
header4 = 'Time took';
fprintf(fid, [header1  '    ' header2 '      ' header3 '        ' header4 '\n']);
% same for adjacency
communities = {}; 
for k = 3 : NUMBER_K_COMMUNITIES; 
  % Use kmeans to cluster the points in the vector W  
  [classes, centers, sumd, D, timeElapsed] = kmeans(W, k);
  % Use the columns of A directly for the clustering
  %[classes, centers, sumd, D, timeElapsed] = kmeans(Adj, k);
  
  for kk = 1 : k; 
    %for each community - create a vector with vertices in it
    communities{kk} = nodes(classes == kk, 1);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Compute modularity per each  k kmeans clustering %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Q = 0;
  edge_cuts = 0;
  overallCommunityEdges = 0;
  for i = 1 : length(communities);
    community = communities{i};
    adjSubgraph = Adj(community, community);   
    %get the number of edges
    %divide by 2 as is undirected graph, needed for edge cut calculation
    communityEdges = sum(sum(adjSubgraph))/2;     
    %sum up all edges per partition, needed for edge cut calculation
    overallCommunityEdges = overallCommunityEdges + communityEdges;   
    Q = Q + newman_modularity(Adj, community);
  end
 
  % number of edge cuts per partition NUMBER_EDGES-overallCommunityEdges
  edge_cuts = NUMBER_EDGES-overallCommunityEdges;
  fprintf(fid, '%f             %f              %f            %f \n', [k Q edge_cuts timeElapsed]');
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Find communities using modularity criteria  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partitionModularityCriteria = community_method(Adj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% compute NMI for two methods  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmiKmeansVsDegreeNeighbours = nmi(partitionModularityCriteria, classes);







