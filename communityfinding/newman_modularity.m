%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes Newman's modularity per cluster. In order to compute the Newman's   %
% modularity for the entire partition , all clusetrs resuls needs be added up  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q] = newman_modularity(Adj, comunity)  
    numberEdges = sum(sum(Adj))/2;
    adjSubgraph = Adj(comunity, comunity);
    
    %get the number of edges
    %divide by 2 as is undirected graph
    communityEdges = sum(sum(adjSubgraph))/2;  
   
    %fraction of edges for this community
    edge_fraction = communityEdges / numberEdges;
    
    %expected fraction of internal edges
    edge_fraction_internal = (sum(sum(Adj(comunity, :))) / numberEdges - edge_fraction)^2;
    
    Q = (edge_fraction - edge_fraction_internal);
end