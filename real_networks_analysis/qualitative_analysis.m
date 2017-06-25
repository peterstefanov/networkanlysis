%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Real world network graph  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load edge list of newqtork - textfile in the format 
%Directed edge A->B, where A - FromNodeId	and B - ToNodeId.
EdgeList = load('Wiki-Vote.txt'); 

%set to 1 if dealing with directed graph dataset, otherwise to 0
isDirected = 0;

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
% calculate the degree per vertex and populate the vector
switch (isDirected)
   case 1 %directed graph
     for d = 1 : NUMBER_EDGES; 
      if (~isempty(find(nodes == EdgeList(d, 1))) == 1)
         Adj(find(nodes == EdgeList(d, 1)), find(nodes == EdgeList(d, 2))) = 1; 
      end 
     end 
   otherwise %undericted graph
     for d = 1 : NUMBER_EDGES; 
         Adj(find(nodes == EdgeList(d, 1)), find(nodes == EdgeList(d, 2))) = 1; 
         Adj(find(nodes == EdgeList(d, 2)), find(nodes == EdgeList(d, 1))) = 1;
     end    
endswitch 

%make adjacency matrix sparse
Adj = sparse(Adj);

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
deg = sort(deg, "descend");
deg = full(deg);
for i = 1 : NUMBER_NODES; 
  if (deg(i) >= 50)
    x(i) = log(deg(i));
    y(i) = log(i/NUMBER_NODES);
  end 
end

figure(2);
plot (x, y, 'o');
hold on
title('Wiki Vote log-log plot of the degree rank vs the degree size');
xlabel ("LOG DEGREE");
ylabel ("LOG RANK");

% Count how many data points we have
m = length(x);
% Add a column of all ones (intercept term) to x
X = [ones(m, 1) x(:)];
% Calculate theta
theta = (pinv(X'*X))*X'*y(:);

% Plot the fitted equation from the regression
hold on; % keeps previous plot visible
x_points = X(:,2);
y_points = X*theta;
plot(x_points, y_points, '-')
legend('Degree distribution', 'Linear regression')
hold off 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Calculate the slope of line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second parameter return for theta contains the slope value, still no harm to
%calculate the slope manually
y1 = y_points(2);
y2 = y_points(m - 1);

x1 = x_points(2);
x2 = x_points(m - 1);

slope = (y2 - y1) / (x2 -x1);

alpha = abs(slope) + 1;