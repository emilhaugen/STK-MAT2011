function lasso_pred = lasso_predict(target, D, L, lambda, BLOCK_LEN)
    % Use a trained dictionary to predict target.
    %   Lasso regression is run on each block of the target image to
    %   recosntruct it from the atoms of the dictionary D.
    %   Assuming target and blocks are square.
    %
    %   param target (matrix): image to predict patch-wise
    %
    %   param D (matrix): dictionary
    %
    %   param L (matrix): sigma^(-.5)
    %
    %   param lamdba (float): lasso regularization parameter
    %
    %   param BLOCK_LEN (int): length of square blocks 
    %
    
    
    SUBSET_LEN = length(target(1,:));
    lasso_pred = zeros(SUBSET_LEN, SUBSET_LEN);
    N = SUBSET_LEN/BLOCK_LEN;
    D_new = L*D; %  decorrelate
    tic
    for t = 0:(N-1)
        disp(t);
        for k = 0:(N-1)
            rows = t*BLOCK_LEN+1:(t+1)*BLOCK_LEN;
            columns = k*BLOCK_LEN+1:(k+1)*BLOCK_LEN;
            patch = target(rows, columns);
            beta = lasso(D_new, L*patch(:), "Lambda", lambda);
            reconstruct_patch = reshape(D_new*beta, BLOCK_LEN, BLOCK_LEN);
            lasso_pred(rows, columns) = reconstruct_patch;
        end
    end  
    toc
end

