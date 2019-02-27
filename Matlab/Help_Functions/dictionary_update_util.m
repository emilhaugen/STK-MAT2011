function [D, D_prev, updated] = dictionary_update_util(D_prev, A, B)
    % DICTIONARY_UPDATE_UTIL utility function called by update_dictionary()
    %   
    %  See dictionary_update() for detaild explanation of input args. 
    %
    %  Return D (matrix): updated dictionary matrix
    %
    %  Return D_prev (matrix): input dictionary matrix, unchanged.
    %
    %  Return updated (1D array): vector of length equal to no. of 
    %    columns in dictionary matrix. Has 1 in position j if column
    %    j was updated, 0 otherwise. (Mainly for testing/debugging).
    %
    
    %{
    diag(A))) < 1e-10
    error("Error: Zero diagonal element in A")
    %}
    CODE_LEN = length(D_prev(1,:));
    D = D_prev;
    updated = zeros(CODE_LEN, 1);
    
    for j = 1:CODE_LEN %iterate over columns of D
        %update column j in D iff A(j,j) != 0.
        if abs(A(j,j)) > 1e-8
            updated(j) = 1;
            dj = B(:,j) - D_prev*A(:,j) + D_prev(:,j)*A(j,j); 
            D(:,j) = dj / (norm(dj) * A(j,j));
        end
    end 
end

