function [D, D_prev, updated] = dictionary_update_util(D_prev, A, B)
    % DICTIONARY_UPDATE_UTIL utility function called by update_dictionary()
    %   
    %  param D_prev (matrix): current version of dictionary matrix
    %  
    %  param A (matrix): outer product X * X' from most recent LASSO
    %
    %  param B (matrix): outer product U * X' from most recent LASSO
    %
    %   
    %  Perform one iteration of block coordinate descent, i.e., loop
    %   over all columns of dictionary ONCE. Only update column j of D
    %   when A(j,j) is non-zero, to avoid zero division.
    %
    %
    %  return D (matrix): updated dictionary matrix
    %
    %  return D_prev (matrix): input dictionary matrix, unchanged. To check
    %                            convergence in dictionary_update()
    %
    %  return updated (1D array): vector of length equal to no. of 
    %                               columns in dictionary matrix. 
    %                               Has 1 in position j if column
    %                               j was updated,   0 otherwise. 
    %                               Mainly for testing/debugging.
    %
    
    CODE_LEN = length(D_prev(1,:)); % no. of columns in dictionary 
    D = D_prev;
    updated = zeros(CODE_LEN, 1);
    
    for j = 1:CODE_LEN % iterate over columns of D
        % update column j in D iff A(j,j) != 0.
        if abs(A(j,j)) > 1e-8
            updated(j) = 1;
            dj = B(:,j) - D_prev*A(:,j) + D_prev(:,j)*A(j,j); 
            D(:,j) = dj / (norm(dj) * A(j,j));
        end
    end 
end
