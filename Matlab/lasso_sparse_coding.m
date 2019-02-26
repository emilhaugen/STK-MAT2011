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
    if abs(min(diag(A))) < 1e-10
        error("Error: Zero diagonal element in A")
    
    
    for j = 1:T
        [b, fitinfo] = lasso(D, U(:,j)); 
        %b is matrix, columns correspond to distinct lambda.
        [minMSE, minIndex] = min(fitinfo.MSE); %use lambda which min. MSE
        %fitinfo.MSE(minIndex)
        X(:,j) = b(:,1);
    end 
end
