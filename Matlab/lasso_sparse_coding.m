function X = lasso_sparse_coding(U, D, X, T)
    %Get sparse codes x^t for vectors u^t (columns in U) by solving 
    %LASSO with least angle regression (LARS), keeping D fixed,
    %
    
    for j = 1:T
        [b, fitinfo] = lasso(D, U(:,j)); 
        %b is matrix, columns correspond to distinct lambda.
        [minMSE, index] = min(fitinfo.MSE); %use lambda to min. MSE
        fitinfo.MSE(index); 
        X(:,j) = b(:,index);
    end 
end
