function community = community_method(adj)    
    tic;
    % Number of edges 
    numberEdges = sum(sum(adj))/2;   
    % Degree vector
    degree = sum(adj,2);
    % Community degree
    community_degree = degree;    
    % Community per node
    community = 1 : length(adj);

    % Neighbours list for each node
    for i = 1 : length(adj)
        Nbs{i} = find(adj(i, :));
        Nbs{i}(Nbs{i} == i) = [];
    end
       
    % Until nodes and communities can be updated
    nodes_flag = true;
    communities_flag = true;
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
                idx = randi(length(list_of_nodes));
                node = list_of_nodes(idx);
                list_of_nodes(idx) = [];
                % Find neighbour communities of node
                node_n_community = unique(community(Nbs{node}));
                node_n_community(node_n_community == community(node)) = [];
                % For each neighbour of community of node n
                best_criteria = 0;
                neighbours_of_node = Nbs{node};
    
                nb1 = neighbours_of_node(community(neighbours_of_node) == community(node));
           
                sum_nb1 = -sum(adj(node, nb1));
                community_degree_node = community_degree(community(node)) - degree(node);
                for i = 1 : length(node_n_community)
                    current = node_n_community(i);
                    next_neighbour_of_node = neighbours_of_node(community(neighbours_of_node) == current);
                    criteria = sum_nb1 + sum(adj(node, next_neighbour_of_node));
                    criteria = (criteria + (degree(node)*(community_degree_node - community_degree(current)))/(numberEdges * 2))/numberEdges;
                    if criteria > best_criteria
                        best_criteria = criteria;
                        best_node_community = node_n_community(i);
                    end
                end

                if best_criteria > 0
                    % Update total degree of communities
                    community_degree(community(node)) = community_degree(community(node)) - degree(node);
                    community_degree(best_node_community) = community_degree(best_node_community) + degree(node);
                    % Update community of the node with the best possible
                    community(node) = best_node_community;
                    % set flags to true
                    merge_flag = true;
                    communities_flag = true;
                end

            end

        end 
        nodes_flag = false;
            
        if ~communities_flag
            break;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Merging communiites %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Community merging
        merge_flag = true;
        while merge_flag
            %fprintf('Community loop\n');
            merge_flag = false;

            % Create community list community_list
            community_list = unique(community);

            % While the list of candidates is not finished
            while ~isempty(community_list)

                % Get a random community n from the community_list and remove it from it
                idx = randi(length(community_list));
                community_n = community_list(idx);
                community_list(idx) = []; 

                % Find neighbour communities of community_n
                neighbour_community_n = find(community == community_n);
                node_n_community = unique(community(unique([Nbs{neighbour_community_n}])));
                node_n_community(node_n_community == community_n) = [];
                    
                % For each neighbour community of community_n
                best_criteria = 0;
                sum_dn1 = sum(degree(neighbour_community_n));
                for ncom_idx = 1 : length(node_n_community)
                    % Compute criteria for merging community_n with current community
                    next_neighbour_community_n = community == node_n_community(ncom_idx);
                    criteria = (sum(sum(adj(neighbour_community_n, next_neighbour_community_n))) - sum_dn1*sum(degree(next_neighbour_community_n))/(numberEdges * 2))/numberEdges;
                    if criteria > best_criteria
                        best_criteria = criteria;
                        best_node_n_community = node_n_community(ncom_idx);
                    end
                end
                % If a move is worth it, do it
                if best_criteria > 0
                    % Update total degree of communities
                    community_degree(best_node_n_community) = community_degree(best_node_n_community) + community_degree(community_n);
                    community_degree(community_n) = 0;
                    % Merge communities
                    community(neighbour_community_n) = best_node_n_community;
                    % set flags to true
                    merge_flag = true;
                    nodes_flag = true;
                end
            end
        end 
        communities_flag = false;
    end 
    
    % Update/rearange communities
    temporary_community = unique(community);
    for i=1:length(community)
        community(i) = find(temporary_community == community(i));
    end
    community = community';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% compute Newman's modularity and edge cut %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen('newmans.txt', 'w');
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
    % number of edge cuts per partition NUMBER_EDGES-overallCommunityEdges
    edge_cuts = numberEdges-overallCommunityEdges;
    fprintf(fid, '%f             %f              %f            %f \n', [numberCommunities modularity edge_cuts timeElapsed]');
    fclose(fid);
end  