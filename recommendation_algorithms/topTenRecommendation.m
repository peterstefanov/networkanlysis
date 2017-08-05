Y = load('rating.txt');
Y = [Y(:, 1:2), Y(:,4)];

% get all unique user, sorted
users = sort(unique([Y(:, 1)]));
% number of users 
NUMBER_USERS = numel(users);

fid = fopen('resultsWithoutTrustDataFirstTenUsers.txt', 'w');
header1 = 'User Id ';
header2 = 'ProductId';
fprintf(fid, [header1  '      ' header2 '      ' header2 '      ' header2 '     ' header2 '       ' header2 '       ' header2 '        ' header2 '       ' header2 '       ' header2 '       ' header2 '\n']);

%sparse metrix representation users - rows  items column
trainSet = sparse(Y(:,1), Y(:,2), Y(:,3));
tic;     
for activeUser = 1 : NUMBER_USERS, %10 - for ten users only 
      
    %get all products rated by the activeUser       
    %find indices to elements in user column
    indActiveUser = Y(:, 1) == activeUser;
    %get the submatrix with the active user only
    submatrixActiveUser = Y(indActiveUser, :);
    %vector with products active user rated
    activeUserItems = submatrixActiveUser(:, 2);
    
    %number of neighbours(most similar to the active user) to be considered
    k = 10;
              
    %similarities between active user and all the rest,
    %based on rating the same/most items like active user
    [sim, user_row] = computeSimilarities(activeUser, trainSet);
    
    %eliminate the active user
    sim(activeUser) = -1;          
        
    %sort similarities in descending order, means the users rated
    %most items adn with highest rates will be on top      
    [s, indx] = sort(sim, 'descend');

    %set of users id most similar to me
    neighbours = indx(1: k);  
        
    %vector holds all neighbours top ten rated products
    neighboursTopRatedProducts = []; 
        
    %for all k neighbours get the top 10 rated products 
    for i = 1 : length(neighbours),
      nextNeighbour = neighbours(i);
      %get only high rated products, rate > 4 or product rated multiple times
      allRatedProducts = trainSet(nextNeighbour, :);
      maskOnlyHighRated = trainSet(nextNeighbour, :) > 4;
      allRatedProducts(maskOnlyHighRated == 0) = -1;
      
      [nextNeighbourProducts, itemsIndex] = sort(allRatedProducts, 'descend');
      itemsIndex = reshape(itemsIndex, numel(itemsIndex), 1);
      nextNeighbourTopRatedProducts = itemsIndex(1: 10);
      %add next neighbour ten top rated products to the vector
      neighboursTopRatedProducts = cat(1, neighboursTopRatedProducts, nextNeighbourTopRatedProducts);         
    end
        
    %count how many times each product appears
    x = unique(neighboursTopRatedProducts);
    N = numel(x);
    count = zeros(N, 2);
    for k = 1 : N,
      %ignore products user already rated
      if any(activeUserItems == x(k)) > 0,
         count(k, 1) = x(k);
         count(k, 2) = -1;
      else   
         count(k, 1) = x(k);
         count(k, 2) = sum(neighboursTopRatedProducts == x(k));
      end   
    end
        
    offset = 9;     
    if N < 20,
      offset = 1;
    end  
        
    %get best top ten rated products from all neighbours
    topRatedProducts = sortrows(count, 2);    
    topRatedProducts = flipud(topRatedProducts((N - offset):end,:));
    top = topRatedProducts(:, 1);  
        
    fprintf(fid, '%f      %f     %f     %f    %f      %f      %f      %f      %f      %f      %f      \n', [activeUser top(1) top(2) top(3) top(4) top(5) top(6) top(7) top(8) top(9) top(10)]');
end  
fclose(fid);
timeElapsed = toc;