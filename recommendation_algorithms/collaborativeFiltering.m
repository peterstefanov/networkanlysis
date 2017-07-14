
tic;
%Y = load('ml.dat');

Y = load('rating.txt');
Y = [Y(:, 1:2), Y(:,4)];

% get all unique user, sorted
users = sort(unique([Y(:, 1)]));
% number of users 
NUMBER_USERS = numel(users);

% get all unique user, sorted
products = sort(unique([Y(:, 2)]));
% number of users 
NUMBER_PRODUCTS = numel(products);

%rust = load('trust.txt');
%trust = [Y(:, 1:2), Y(:,4)];

fid = fopen('C:\Users\stefp\COMP42270\networkanalysisgit\networkanlysis\recommendation_algorithms\allPerUser.txt', 'w');
header1 = 'User Id ';
header2 = 'Product Id';
header3 = 'Avg ratting';
fprintf(fid, [header1  '             ' header2 '            ' header3 '\n']);

% random permutation 
%p = randperm(length(Y));

%reorganzing 
%Y(:,1) = Y(p,1);
%Y(:,2) = Y(p,2);
%Y(:,3) = Y(p,3);

% split into 5 sets
numTrans = length(Y);

first=1;
testSize = floor(numTrans/200);
last = testSize;
sumAbsErr = 0;
totalTrans = 0;
%matrix to store the avarage ratting for product per user
usersNeighboursRattingPerProduct = sparse(NUMBER_USERS, NUMBER_PRODUCTS);
%for i=1:200,
    
    testY = Y(first:last,:);
    trainY = [Y(1:(first-1),:);Y((last+1):end,:)];
    
    %first = first+testSize;
    %last = last+testSize;

    %sparse metrix representation users - rows  items column
    %trainSet = sparse(trainY(:,1),trainY(:,2),trainY(:,3));
    trainSet = sparse(Y(:,1),Y(:,2),Y(:,3));

    testSet = sparse(testY(:,1),testY(:,2),testY(:,3));
    
    
    %for trans=1:length(testY),
    for trans=1:length(testY),  
        activeUser = testY(trans,1);
        activeItem = testY(trans,2);
        activeRating = testY(trans,2);
        
        k = 20;
        
        %how many items user rated
        numberItemsRatedByActiveUser = numel(find(activeUser == Y(:, 1)));       
        %matrix to store top k users that rated each item like me,
        %needed to decide which users are more simmilar to me,
        %based on users rated most items like me
        itemsUsersRatedLikeMe = sparse(numberItemsRatedByActiveUser, k);
             
        user_row = trainSet(activeUser,:);     
        %similarities between active user and all the rest,
        %based on if they rate the same item like active user
        sim = computeSimilarities(activeUser, trainSet);
        
        %sim = 0.5*sim +0.5*trustToUser
        %eliminate myself
        sim(activeUser) = -1;
                      
        %binary variable false or true 
        mask = trainSet(:, activeItem) > 0;
        %eliminate not rated similarities
        sim(mask == 0) = -1;  
        
        [s, indx] = sort(sim, 'descend');
        break;
        %set of users id most similar to me
        neighbours = indx(1:k);

        mask=sim(neighbours)>0;
              
        neighbours = neighbours(mask);
        
        if (length(neighbours) == 0)
            predictRating = 0;
        else

        %get all users rated like me the active item 
        ratedLikeMe = find(activeRating == trainSet(neighbours, activeItem));
  
        break;
        %avarage reating that the neighbours giving for this product/item
         predictRating = round(mean(trainSet(neighbours,activeItem)));

        end
        rattingCompareWithNeighbours  = predictRating - activeRating;
        usersNeighboursRattingPerProduct(activeUser, activeItem) = predictRating; 
        
       
        %how well the algorithm has done
       sumAbsErr = sumAbsErr+abs(predictRating - activeRating); 
       totalTrans = totalTrans+1;
       
       %if (mod(totalTrans,100)==0)          
           %fprintf('%e\n', sumAbsErr/totalTrans);
       %end   
  fprintf(fid, '%f             %f              %f \n', [activeUser activeItem rattingCompareWithNeighbours]');
  
  
  end  
   fclose(fid);
%end

MAE = sumAbsErr/totalTrans;
% record time taken
timeElapsed = toc;