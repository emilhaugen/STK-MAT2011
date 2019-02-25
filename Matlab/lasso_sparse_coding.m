function X = lasso_sparse_coding(U, D)
    %Get sparse codes for vectors u^t (columns in U) by solving 
    %LASSO with least angle regression (LARS), keeping D fixed
    
    Udim = size(U);
    k = Udim(1);
    for j = 1:2
        b = lasso(D, U(:,j));
    end 
    X = b;
end
