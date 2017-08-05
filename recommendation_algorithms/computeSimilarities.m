function [sim, user_row] = computeSimilarities(user, trainSet)

user_row = trainSet(user,:);
%Active user's rating multiply by the ratings other users gave for the same 
%products active user rated. Consider: if active user rated product 1 with 5 
%and some user rated with 3, it will add 15 for this product, if any of the 
%other products that the active user rated apperas to berated from the same 
%user the same will be done for it and will be added to 15. Example user one 
%rated products 1 to 10 with (2 2 2 5 3 4 1 4 4 1), and user 1174 rated four 
%of his products with 1 with 3; 7 with 4; 9 with 4; 8 with 5, so the overall 
%result will be 3*2 + 4*1 + 4*4 + 5*4 = 46, therefore the user with the highest 
%number is more similar to the active user. This user rated four of active user 
%products and have score of 46.
sim = trainSet * user_row';