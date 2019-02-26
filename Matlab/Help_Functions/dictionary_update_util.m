function [D, D_prev] = dictionary_update_util(D_prev, A, B)
    % DICTIONARY_UPDATE_UTIL utility function called by update_dictionary()
    %   
    % 
    % Return D (matrix): updated dictionary matrix
    %
    % Return D_prev (matrix): input dictionary matrix, unchanged.
    %
    
    if abs(min(diag(A))) < 1e-10
        error("Error: Zero diagonal element in A")
    
    CODE_LEN = length(D_prev(1,:));
    D = D_prev;
    for j = 1:CODE_LEN %iterate over columns of D
        dj = B(:,j) - D_prev*A(:,j) + D_prev(:,j)*A(j,j); 
        D(:,j) = dj / (norm(dj) * A(j,j));
    end 
end

