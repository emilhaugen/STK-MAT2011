%
%  Script to test dictionary_update function. Not much to see.
%  
%   Manually testing alternating optimization between 
%       LASSO and BCD (block coord. descent)
%
%
%
addpath("Help_Functions"); %make help functions available
addpath("Preprocess");

%set constants
FILENAME = "Restoration.txt";
SUBSET_LEN = 64;
BLOCK_LEN = 8; %data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; %length of sparse code vectors, no. of colums in dictionary

%read ascii data
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);
Udim = size(U);
T = Udim(2); %using notation T for no. of data vectors as in paper

%randomly initialize dictionary
D = init_dict(U, CODE_LEN);

%get initial LASSO estimates for X
X = zeros(CODE_LEN, T);
X = lasso_sparse_coding(U, D, X, T);
A = sum_of_outer_prods(X, X);
B = sum_of_outer_prods(U, X);

%test update_deictionary() function
tol = 1e-9;
max_iter = 50;
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, max_iter);


%evaluate results of dictionary update  
histogram(D)
max(max(D))
norm(D, "fro")
norm(D - D_prev, "fro")
nnz(updated)
num_iter

%repeat LASSO
X = lasso_sparse_coding(U, D, X, T);
A = sum_of_outer_prods(X, X);
B = sum_of_outer_prods(U, X);

%repeat update_dictionary()
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, max_iter);

%evaluate results of dictionary update  
histogram(D)
max(max(D))
norm(D, "fro")
norm(D - D_prev, "fro")
nnz(updated)
num_iter

%{
norm(D - D_prev, "fro")
max(max(D))
max(max(D_prev))
z = zeros(T);
z = z(:, 1);
for p = 1:T
    z(p) = CODE_LEN - nnz(X(:,p));
end    
z    
%}




