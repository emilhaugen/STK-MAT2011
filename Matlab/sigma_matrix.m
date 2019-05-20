function [sigma] = sigma_matrix(C, M)
    %SPATIAL_COVARIANCE get spatial covariance 
    %  extract spatial covariance
    %   param C (matrix): spatial covariance matrix
    %   param M (int): side length of square blocks
    
    N = length(C(1,:)); % C is square
    sigma = zeros(M*M);
    
    % temp. variable to manage indices
    ind = zeros(M*M, 2);
    for i = 0:(M-1)
        ind(M*i + 1 : M*(i+1), 1) = i+1;
        ind(M*i + 1 : M*(i+1), 2) = 1:M;
    end   
    for i = 1:M*M
        for j = 1:M*M
            row_ind = ind(i,1) - ind(j,1) + 1;
            col_ind = ind(i,2) - ind(j,2) + 1; 
            if row_ind < 1 % circular matrix
                row_ind = N + row_ind;
            end
            if col_ind < 1
                col_ind = N + col_ind;
            end 
            sigma(i,j) = C(row_ind, col_ind); % col row mixed up
        end
    end
end

