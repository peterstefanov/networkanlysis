%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Price's %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NUMBER_NODES = 1024;
%Adj = priceModel(NUMBER_NODES);

%%make adjacency matrix sparse
%Adj = sparse(Adj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Kronecker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCALE = 10;
edgefactor = 16;

kron = kroneckerModel(SCALE, edgefactor);
Adj = sparse (kron(1,:)+1, kron(2,:)+1, ones (1, size (kron, 2)));
% remove self-loops
Adj = Adj - diag(diag(Adj));
NUMBER_NODES = length(Adj);
%%%%%%%%%%%%%%%%%%%%% END GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no multiple edges
Adj = Adj > 0; 

indeg = sum(Adj);
outdeg = sum(Adj');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Compute clustering coefficient %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multiplication change in the clustering coefficient formula
% to reflect the edge count for directed/undirected graphs
if Adj == transpose(Adj); %undirected 
   deg = indeg + diag(Adj)';  
   coeff = 2;
else %directed
   deg = indeg + outdeg; 
   coeff = 1; 
end

% initialize local clustering coefficient vector
Ci = zeros(NUMBER_NODES, 1); 

for i = 1 : NUMBER_NODES  
  if deg(i) == 1 || deg(i) == 0; 
    Ci(i) = 0; 
    continue; 
  end

  neighbors = find(Adj(i, :) > 0); 
  subgraph = Adj(neighbors, neighbors);
    
  if (subgraph == transpose(subgraph));   
    edges = sum(sum(subgraph)) / 2;    
  else   
    edges = sum(sum(subgraph));   
  end

  Ci(i) = coeff * edges / (deg(i) * (deg(i) - 1));
end

% average clustering coefficient
C = sum(Ci) / NUMBER_NODES;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Compute Assortativity  coefficient %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute degrees
Adj = double(Adj ~= 0);

if Adj == transpose(Adj); %undirected 
   alldeg = sum(Adj);
   [i, j] = find(triu(Adj, 1) > 0);   
else %directed 
   % indegree - column sum of Adj
   indegree = sum(Adj, 1);  
   % outdegree - row sum of Adj
   outdegree = sum(Adj, 2)'; 
   % degree = indegree + outdegree  
   alldeg = indegree + outdegree;     
   [i, j] = find(Adj > 0);
end
   K = length(i);

for k = 1 : K
    degi(k) = alldeg(i(k));
    degj(k) = alldeg(j(k));
end;

% compute assortativity
r = (sum(degi.*degj) / K - (sum(0.5*(degi + degj)) / K)^2) / (sum(0.5*(degi.^2 + degj.^2)) / K - (sum(0.5*(degi + degj)) / K)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Compute Degree Distribution  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
spy (Adj);
hold on
title('Synthetic graph - density plot');
deg = sort(deg, "descend");

for i = 1 : NUMBER_NODES; 
  if (deg(i) >= 50)
    degree_plot(i) = log(deg(i));
  end 
end