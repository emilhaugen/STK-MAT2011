function data = patches_to_original(patches, BLOCK_LEN, ncol, nrow)
    % PATCHES_TO_ORIGINAL convert data on patch form back to original
    %                       matrix form. Inverse of ascii_to_data_matrix().
    %
    %  param patches (matrix): data with patches stored in columns
    %
    %  param BLOCK_LEN (matrix): the side lenght of each square patch 
    %
    %  param ncol (int): no. of columns in original format
    %
    %  param nrow (int) no. of rows in original format
                        
    if mod(ncol, BLOCK_LEN) ~= 0 || mod(nrow, BLOCK_LEN) ~= 0
        error("BLOCK_LEN must divide both ncol and nrow");
    end
    data = zeros(nrow, ncol);
    
    counter = 1; % count which block is being processed
    for i = 1:(nrow/BLOCK_LEN)
        row_start = 1 + (i-1)*BLOCK_LEN;
        row_end = i*BLOCK_LEN;
        for j = 1:(ncol/BLOCK_LEN)
             col_start = 1 + (j-1)*BLOCK_LEN;
             col_end = j*BLOCK_LEN;
             data(row_start:row_end,col_start:col_end) = ...
              transpose(reshape(patches(:,counter), BLOCK_LEN, BLOCK_LEN));
             counter = counter + 1;
        end
    end    
end

