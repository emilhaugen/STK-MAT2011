function data = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN)
    % ASCII_TO_DATA_MATRIX convert ASCII file to data matrix U used in DL.
    %   
    %  pram FILENAME: (string) containing location of ASCII file
    %  
    %  param SUBSET_LEN: (int) usually a power of two
    %
    %  param BLOCK_LEN: (int) must divide SUBSET_LEN 
    % 
    %    Divides subset (square) image into blocks of size 
    %        BLOCK_LEN*BLOCK_LEN 
    %       
    %    Each block is saved as a column of a matrix data.
    %     
    
    U = dlmread(FILENAME);
    U_subset = U(1:SUBSET_LEN, 1:SUBSET_LEN); %top left corner of image
    %imshow(mat2gray(U_subset))

    subset_block_ratio = SUBSET_LEN / BLOCK_LEN; %no. of blocks per column/row
    p = subset_block_ratio^2; %no. of columns in final data matrix

    M = BLOCK_LEN^2; %dimension of data vectors in DL
    data = zeros(M, p); %initialize data matrix being sent to DL
    counter = 1; %count which block is being processed

    for i = 1:subset_block_ratio
        row_start = 1 + (i-1)*BLOCK_LEN;
        row_end = i*BLOCK_LEN;
        for j = 1:subset_block_ratio
             col_start = 1 + (j-1)*BLOCK_LEN;
             col_end = j*BLOCK_LEN;
             block = U_subset(row_start:row_end,col_start:col_end);
             data(:,counter) = flatten(block);
             counter = counter + 1;
        end
    end
end

