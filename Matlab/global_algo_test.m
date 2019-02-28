%
%  Global script to test convergence.
%
%set constants
addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 64; %should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 8; %data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; %length of sparse code vectors, no. of colums in dictionary
tol = 1e-6;
max_iter = 5;
dict_iter = 50;

%read ascii data
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);
%randomly initialize dictionary
D = init_dict(U, CODE_LEN, false);

X = lasso_sparse(U, D);

