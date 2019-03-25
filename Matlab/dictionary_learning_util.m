function [D, D_prev, X, A, B, updated, i] = dictionary_learning_util(U, D_prev, lambda, dict_tol, dict_iter)
    % DICTIONARY_LEARNING_UTIL utility function called by
    %                           dictionary_learning()
    %
    %   Perform LASSO regression using lasso_sparse_coding(), 
    %        then do dictionary update with dictionary_update().
    %
    %
    %   param U (matrix): input data to DL algorithm from ascii file
    %
    %   param D_prev (matrix): current version of dictionary matrix.
    %
    %   param lambda (float): fixed parameter passed to LASSO
    %
    %   return D (matrix): dictionary learned after one LASSO and
    %                        one dictionary update iteration
    %
    %   return D_prev (matrix): unchanged input dictionary   
    %
    %   return X (matrix): sparse code vectors after one LASSO iteration
    %   
    %   return A (matrix): X * X'
    %
    %   return B (matrix): U * X' 
    
    
    % perform LASSO
    [X, A, B] = lasso_sparse_coding(U, D_prev, lambda);
    
    % iterative dictionary update subroutine
    [D, D_prev, updated, i] = dictionary_update(D_prev, A, B, dict_tol, dict_iter);
    %norm(D - D_prev, "fro");

end

