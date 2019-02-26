function X = lasso_sparse_coding(U, D, X, T)
    %LASSO_SPARSE_CODING Get sparse codes x^t for vectors u^t 
    %                     (columns in U) by solving LASSO, keeping D fixed.
    %
    %  param U (matrix): data matrix from ascii files
    % 
    %  param D (matrix): current dictionary matrix 
    %
    %  param X (matrix): current sparse code vectors
    %
    %  param T (int): no. of columns in U, aka no. of data vectors
    
    for j = 1:T
        %b is matrix, columns correspond to distinct lambda.
        [b, fitinfo] = lasso(D, U(:,j)); 
        [minMSE, minIndex] = min(fitinfo.MSE); %use lambda which min. MSE
        %fitinfo.MSE(minIndex)
        X(:,j) = b(:,minIndex);
    end 
end
