function X = lasso_sparse_coding(U, D, X, T)    
    % LASSO_SPARSE_CODING Get sparse codes x^t for vectors u^t 
    %                     (columns in U) by solving LASSO, keeping D fixed.
    %
    %  param U (matrix): data matrix from ascii files
    % 
    %  param D (matrix): current dictionary matrix 
    %
    %  param X (matrix): current sparse code vectors. All of X is reset
    %                      in each iteration of LASSO, X is input here
    %                      only to avoid repeating dimension calculation.
    %
    %  param T (int): no. of columns in U, aka no. of data vectors
    %
    %  Solve LASSO regression problem with constant matrix D and response 
    %   U(:,j) and save resulting coefficients as X(:,j) for j = 1,...,T.
    %   Uses coefficients corresponding to minimum MSE.
    %   
    %  return X (matrix): with LASSO coefficients as columns.
    
    for j = 1:T
        %b is matrix, columns correspond to distinct lambda.
        [b, fitinfo] = lasso(D, U(:,j)); 
        [minMSE, minIndex] = min(fitinfo.MSE); %use lambda which min. MSE
        %fitinfo.MSE(minIndex)
        X(:,j) = b(:,minIndex);
    end 
end
