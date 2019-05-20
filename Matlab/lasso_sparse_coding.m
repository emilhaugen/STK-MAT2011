function [X, A, B] = lasso_sparse_coding(U, D, lambda)    
    % LASSO_SPARSE_CODING Get sparse codes x^t for vectors u^t 
    %                     (columns in U) by solving LASSO, keeping D fixed.
    %
    %  param U (matrix): data matrix from ascii files
    % 
    %  param D (matrix): current dictionary matrix 
    %
    %  Solve LASSO regression problem with constant matrix D and response 
    %   U(:,j) and save resulting coefficients as X(:,j) for j = 1,...,T.
    %   
    %  return X (matrix): with LASSO coefficients as columns.
    
    T = length(U(1,:)); % no. of columns in U, i.e. no. of data vectors
    CODE_LEN = length(D(1,:)); % no. of columns in dictionary
    X = zeros(CODE_LEN, T);
    for j = 1:T
        if mod(j, 50) == 0
            fprintf("Lasso Iteration %d of %d\n", j, T)
        end    
        % b is matrix, columns correspond to distinct lambda
        % or just 1D-vector if lambda specified
        X(:,j) = lasso(D, U(:,j), "Lambda", lambda); 
    end
    A = X*X';
    B = U*X';
end
