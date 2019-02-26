function outer = sum_of_outer_prods(X, Y)
    % SUM_OF_OUTER_PRODS help function to calculate sum of
    %                      outer products of columns in X and Y
    %   
    %   param X (matrix): c x T matrix
    % 
    %   param Y (matrix): m x T matrix
    %
    %   Calculate sum over j=1,..., T of the outer product matrices
    %     X(:,j) * transpose(Y(:,j))                  
    %
    %   Return outer (matrix): Sum of outer products. (c x m matrix)
    %
    outer = X * transpose(Y);
end

