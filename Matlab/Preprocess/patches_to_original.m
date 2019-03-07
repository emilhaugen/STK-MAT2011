function data = patches_to_original(patches)
    % PATCHES_TO_ORIGINAL convert data on patch form back to original
    %                       matrix form. Inverse of ascii_to_data_matrix().
    
    
    [m, T] = size(patches); % patch side length is sqrt(m)
                        % by def, T = (SUBSET_LEN / BLOCK_LEN)^2
                        
    BLOCK_LEN = sqrt(m);                        
    SUBSET_LEN = sqrt(T) * BLOCK_LEN;
    data = zeros(SUBSET_LEN, SUBSET_LEN);

    subset_block_ratio = SUBSET_LEN / BLOCK_LEN;

    counter = 1; %count which block is being processed
    for i = 1:subset_block_ratio
        row_start = 1 + (i-1)*BLOCK_LEN;
        row_end = i*BLOCK_LEN;
        for j = 1:subset_block_ratio
             col_start = 1 + (j-1)*BLOCK_LEN;
             col_end = j*BLOCK_LEN;
             data(row_start:row_end,col_start:col_end) = ...
              transpose(reshape(patches(:,counter), BLOCK_LEN, BLOCK_LEN));
             counter = counter + 1;
        end
    end    
end

