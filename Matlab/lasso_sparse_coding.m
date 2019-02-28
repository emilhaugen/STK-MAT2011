function [X, A, B] = lasso_sparse_coding(U, D)    
    % LASSO_SPARSE_CODING Get sparse codes x^t for vectors u^t 
    %                     (columns in U) by solving LASSO, keeping D fixed.
    %
    %  param U (matrix): data matrix from ascii files
    % 
    %  param D (matrix): current dictionary matrix 
    %
    %  Solve LASSO regression problem with constant matrix D and response 
    %   U(:,j) and save resulting coefficients as X(:,j) for j = 1,...,T.
    %   Uses coefficients corresponding to minimum MSE.
    %   
    %  return X (matrix): with LASSO coefficients as columns.
    
    T = length(U(1,:)); % no. of columns in U, i.e. no. of data vectors
    CODE_LEN = length(D(1,:)); % no. of columns in dictionary
    X = zeros(CODE_LEN, T);
    
    for j = 1:T
        %b is matrix, columns correspond to distinct lambda.
        [b, fitinfo] = lasso(D, U(:,j)); 
        [minMSE, minIndex] = min(fitinfo.MSE); %use lambda which min. MSE
        %fitinfo.MSE(minIndex)
        X(:,j) = b(:,minIndex);
    end 
    A = X*X';
    B = U*X';
end
