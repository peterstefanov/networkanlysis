function [community, timeElapsed] = community_method(adj)    
    tic;
    % Number of edges 
    numberEdges = sum(sum(adj))/2;   
    % Degree vector
    degree = sum(adj, 2);
    % Community degree
    community_degree = degree;    
    % Community per node
    community = 1 : length(adj);

    % Neighbours list for each node
    for i = 1 : length(adj)
        Nieghbours{i} = find(adj(i, :));
        Nieghbours{i}(Nieghbours{i} == i) = [];
    end

    % Until nodes and communities can be updated
    nodes_flag = true;
    %communities_flag = true;
    while nodes_flag
        % node can be merge_flag
        merge_flag = true;
        while merge_flag
            merge_flag = false;
            % vector representation of all nodes in the graph
            list_of_nodes = 1 : length(adj);
            % Until list of nodes is not empty
            while ~isempty(list_of_nodes)
                % Get a random node n from the list_of_nodes and remove it from it
                index = randi(length(list_of_nodes));               
                node = list_of_nodes(index);  

                % Get the list of all remaining nodes without the index/choosen node            
                list_of_nodes(index) = [];
      
                % Get all neighbour's nodes of the index/choosen node which are in the community
                node_n_neighbours = unique(community(Nieghbours{node}));                              
                node_n_neighbours(node_n_neighbours == community(node)) = [];

                % For each neighbour of community of node n
                best_criteria = 0;
                node_n_initial_neighbours = Nieghbours{node};  
                % actual neighbours list of community of node n
                neighbours_list_community = node_n_initial_neighbours(community(node_n_initial_neighbours) == community(node));  
                %number of nodes in community
                sum_neighbours_list_community = sum(adj(node, neighbours_list_community));

                %comunity degree without the degree of the node itself
                community_degree_node = community_degree(community(node)) - degree(node);
 
                % For all neighbours of the node 
                for i = 1 : length(node_n_neighbours)
                    current_node = node_n_neighbours(i);
                    % list of neighbours nodes in the community of the current node
                    current_node_neighbours_list = node_n_initial_neighbours(community(node_n_initial_neighbours) == current_node);
                    % positive if the node(choosen) has more nodes in connections with some nodes from the neighbours list 
                    % basically if the node and some of its neighbours has more connections than the number of node neighbours list, set this node as best possible  
                    criteria = sum(adj(node, current_node_neighbours_list)) - sum_neighbours_list_community;
                    if criteria > best_criteria
                        best_criteria = criteria;
                        best_node = node_n_neighbours(i);
                    end
                end

                % once we have the best node 
                if best_criteria > 0
                    % Update total degree of communities
                    community_degree(community(node)) = community_degree(community(node)) - degree(node);
                    community_degree(best_node) = community_degree(best_node) + degree(node);
                    % Update community of the node with the index of the best possible node 
                    % basically marking each node with the node index for which we found  
                    community(node) = best_node;
                    % set flags to true
                    merge_flag = true;
                end
            end
        end 
        nodes_flag = false;           
    end 
    

    % Update/rearange communities
    unique_community_index = unique(community);
    for i = 1 : length(community)
        community(i) = find(unique_community_index == community(i));
    end

    community = community';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% compute Newman's modularity and edge cut %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen('modularityOptim.txt', 'w');
    header1 = 'Community number ';
    header2 = 'Modularity';
    header3 = 'Number of edges cuts';
    header4 = 'Time took';
    fprintf(fid, [header1  '    ' header2 '      ' header3 '        ' header4 '\n']);
    
    nodes = 1 : length(adj);
    nodes = nodes';
    for index = 1 : length(unique(community));       
       %for each community - create a vector with vertices in it
       communitiesNewman{index} = nodes(community == index, 1);
    end
  
    numberCommunities = length(unique(community));
    modularity = 0;
    edge_cuts = 0;
    overallCommunityEdges = 0;
    for b = 1 : numberCommunities;
      comm = communitiesNewman{b};   
      adjSubgraph = adj(comm, comm);   
      %get the number of edges
      %divide by 2 as is undirected graph, needed for edge cut calculation
      communityEdges = sum(sum(adjSubgraph))/2;     
      %sum up all edges per partition, needed for edge cut calculation
      overallCommunityEdges = overallCommunityEdges + communityEdges;   
      modularity = modularity + newman_modularity(adj, comm);
    end
  
    % record time taken
    timeElapsed = toc;
    % number of edge cuts per partition numberEdges-overallCommunityEdges
    edge_cuts = numberEdges - overallCommunityEdges;
    fprintf(fid, '%f             %f              %f            %f \n', [numberCommunities modularity edge_cuts timeElapsed]');
    fclose(fid);
end  