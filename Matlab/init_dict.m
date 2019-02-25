function D = init_dict(U)
    %   INIT_DICT help function to initialize dictionary from data U
    %   
    %   param U (matrix): matrix of data generated by
    %                       ascii_to_data_matrix()
    %   
    %   Randomly sample columns of U to be columns of D
    %
    %   Return D (matrix): Matrix of same dimension as U
    
    
    rng(10, 'twister');
    s = RandStream('mlfg6331_64');
    
    Udim = size(U);
    k = Udim(2); %no. of cols in U
    
    %sample columns with replacement
    D = datasample(s, U, k, 2, 'replace', true); 
end

