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
SUBSET_LEN = 32;
BLOCK_LEN = 8; %data subset will have (SUBSET_LEN / BLOCK_LEN)^2 vectors
CODE_LEN = 100; %length of sparse vectors, no. of colums in dictionary

% set iteration parameters for dictionary learning
tol = 1e-10;
max_iter = 50;
dict_iter = 50;

% read ascii data and initialize dictionary randomly
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);

% dictionary learning
[D, X, A, B, updated, niter_dict, niter_glob] = ...
    dictionary_learning(U, CODE_LEN, tol, max_iter, dict_iter);










