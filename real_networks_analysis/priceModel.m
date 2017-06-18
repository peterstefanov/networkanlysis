%n number of nodes for the graph to be created
function simulate = priceModel(n)

simulate = zeros(n);
simulate(1, 1) = 1; 
vertices = 1;

while vertices < n
    % attach new vertex 
    vertices = vertices + 1;
    simulate(vertices, vertices) = 1; 
    % every time get in-degree 
    in_degree = sum(simulate); 
    
    avg = 0; 
    for k = 1 : vertices
      %this give us the probability
      pk(k) = numel(find(in_degree == k)) / vertices;
      % average in-degree (per vertex), number of new citations per vertex
      avg = avg + pk(k) * k;
    end
     
    for k = 1 : vertices
        % attach new edges with probability (k+1)pk/(avg+1)
        if rand < (k + 1) * pk(k) / (avg + 1); 
          simulate(vertices, k) = simulate(vertices, k) + 1; 
        end
    end
end
% remove self-loops
simulate = simulate - diag(diag(simulate));  






