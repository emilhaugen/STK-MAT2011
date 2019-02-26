function [D, D_prev] = dictionary_update_util(D_prev, A, B)
    % DICTIONARY_UPDATE_UTIL utility function called by update_dictionary()
    %   
    % 
    %
    
    if abs(min(diag(A))) < 1e-10
        error("Error: Zero diagonal element in A")
    
    CODE_LEN = length(D(1,:));
    D_prev = D;
    for j = 1:CODE_LEN %iterate over columns of D
        dj = B(:,j) - D*A(:,j) + D(:,j)*A(j,j); 
        D(:,j) = dj / (norm(dj) * A(j,j));
    end 
end

