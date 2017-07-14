function sim=computeSimilarities(user, trainSet)


user_row = trainSet(user,:);
%active users rating mutiply by the ratings other users gave for the same products
%if active rated product 1 with 5 and some user rated with 3, it will add 15 for this 
%product, if any of the other products that the active user rated the same will be done 
% for it and will be added to 15.
%example user one rated products 1-10 (2 2 2 5 3 4 1 4 4 1), and user 1174 rated four of his products
% with 1-3 7-4 9-4 8-5, so the overall result will be 3*2 + 4*1 + 4*4 + 5*4 = 46, therefore the user
% with the highest nuber is more simmilar to the active user
sim = trainSet*user_row';