function [D, D_input, updated, i]  = dictionary_update(D_prev, A, B, tol, max_iter)
    % DICTIONARY_UPDATE update dictionary using block coordinate descent
    % 
    %  param D (matrix): current dictionary.
    %  
    %  param A (matrix): sum of outer product of current code vectors.
    % 
    %  param B (matrix): sum of outer products of data vectors
    %                        with current code vectors.
    %
    %  param tol (float): tolerance for when D has converged in  
    %                        Frobenius norm.
    %
    %  param max_iter (int): maximum no. of iterations
    %
    %
    %  return D (matrix): updated dictionary matrix.
    %
    %  return D_input (matrix): unchanged input dictionary to check 
    %                            for convergence in dictionary_learning()
    %
    
    D_input = D_prev;
    [D, D_prev, updated] = dictionary_update_util(D_prev, A, B);
    
    i = 0; % continiue updating dictionary until convergence or max iter 
    while norm(D - D_prev, "fro") > tol && i < max_iter
        [D, D_prev, updated] = dictionary_update_util(D, A, B);
        i = i + 1;
    end           
end

