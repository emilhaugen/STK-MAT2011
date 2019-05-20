%
%  
%  Script to test dictionary_learning function.
%
%

% make functions and data available
addpath("Help_Functions"); 
addpath("Preprocess");
addpath("Data");

% set constants
FILENAME = "Restoration.txt";
SUBSET_LEN = 64;
BLOCK_LEN = 8; %data subset will have (SUBSET_LEN / BLOCK_LEN)^2 vectors
CODE_LEN = 75; %length of sparse vectors, no. of colums in dictionary

% set iteration parameters for dictionary learning
tol = 1e-5;
max_iter = 50;
dict_iter = 20;

% read ascii data and initialize dictionary randomly
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);
D_prev = init_dict(U, CODE_LEN);


%{
dictionary_learning_util() testing
D_prev = init_dict(U, CODE_LEN);
j = 0;
while j < max_iter
    [D, D_prev, X, A, B, updated, i] = dictionary_learning_util(U, D_prev, tol, dict_iter);
    nnz(diag(A))
    j = j + 1;
end
diag(A)
i
norm(D - D_prev, "fro")
%}








