function D = dictionary_update(D, A, B, tol, max_iter)
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
    %  Return D (matrix): updated matrix.
    
    Ddim = size(D);
    m = Ddim(1);
    CODE_LEN = Ddim(2);
    
    D_prev = D;
    for j = 1:CODE_LEN %initial update, iterate over columns of D
        dj = B(:,j) - D*A(:,j) + D(:,j)*A(j,j); 
        D(:,j) = dj / (norm(dj) * A(j,j));
    end 
    
    i = 0; %continiue updating dictionary until convergence or max iter 
    while norm(D - D_prev, "fro") > tol && i < max_iter
        D_prev = D;
        for j = 1:CODE_LEN
            dj = B(:,j) - D*A(:,j) + D(:,j)*A(j,j);
            D(:,j) = dj / (norm(dj) * A(j,j));
        end
        i = i + 1;
    end           
end

