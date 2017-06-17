%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Real world network graph  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load edge list of newqtork - textfile in the format 
%Directed edge A->B, where A - FromNodeId	and B - ToNodeId.
EdgeList = load('Wiki-Vote.txt'); 

%sort the edge list
sortrows(EdgeList);

% number of edges in the graph
NUMBER_EDGES = length(EdgeList);

% get all nodes, sorted
nodes = sort(unique([EdgeList(:, 1) EdgeList(:, 2)]));

% number of nodes in the graph
NUMBER_NODES = numel(nodes);

% initialize adjacency matrix nad make it sparse
Adj = zeros(NUMBER_NODES);
Adj = sparse(Adj);
        
% iterate across all edges and populate adjacency matrix with one's
for d = 1 : NUMBER_EDGES;  
   Adj(find(nodes == EdgeList(d, 1)),find(nodes == EdgeList(d, 2))) = 1;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Compute clustering coefficient %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get node degrees
node_degrees = sum(Adj, 2); 
%get the number of triangles for each node
node_triangles = diag(Adj * triu(Adj) * Adj); 
%local clustering coefficient of each node
C_i = zeros(size(node_degrees));
C_i(node_degrees > 1) = 2 * node_triangles(node_degrees > 1) ./ (node_degrees(node_degrees > 1).*(node_degrees(node_degrees > 1) - 1)); 
%average clustering coefficient of the graph
average_clustering = mean(C_i(node_degrees > 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Compute Assortivity  coefficient %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECTED GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute degrees
Adj = double(Adj ~= 0);
% indegree - column sum of Adj
in_degree = sum(Adj, 1);  
% outdegree  - row sum of Adj
out_degree = sum(Adj, 2)'; 
% degree = indegree + outdegree  
all_deg = in_degree + out_degree;  
    
[i,j] = find(Adj > 0);
K = length(i);
for k = 1 : NUMBER_EDGES
    deg_i(k) = all_deg(i(k));
    deg_j(k) = all_deg(j(k));
end;

% compute assortativity
r_directed = (sum(deg_i.*deg_j)/NUMBER_EDGES - (sum(0.5*(deg_i + deg_j))/NUMBER_EDGES)^2)/(sum(0.5*(deg_i.^2 + deg_j.^2))/NUMBER_EDGES - (sum(0.5*(deg_i + deg_j))/NUMBER_EDGES)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNDIRECTED GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alldeg = sum(Adj);
[ii, jj] = find(triu(Adj, 1) > 0);
KK = length(ii);
for k_ = 1 : KK
    degi(k_) = alldeg(ii(k_));
    degj(k_) = alldeg(jj(k_));
end;

% compute assortativity
r_undirected = (sum(degi.*degj) / KK - (sum(0.5*(degi + degj)) / KK)^2) / (sum(0.5*(degi.^2 + degj.^2)) / KK - (sum(0.5*(degi + degj)) / KK)^2);
