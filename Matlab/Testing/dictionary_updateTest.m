%
%  Script to test dictionary_update function. Not much to see here.
%  
%   Manually testing alternating optimization between 
%       LASSO and BCD (block coord. descent)
%
%
%
addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

%set constants
FILENAME = "Restoration.txt";
SUBSET_LEN = 128; %should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 8; %data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; %length of sparse code vectors, no. of colums in dictionary
tol = 1e-6;
max_iter = 5;
dict_iter = 50;
%read ascii data
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);

%randomly initialize dictionary
D = init_dict(U, CODE_LEN);

%get initial LASSO estimates for X
X = lasso_sparse_coding(U, D);
A = X*X';
B = U*X';
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, dict_iter);
norm(D*X - U, "fro")
num_iter
nnz(updated)
D1 = D;

X = lasso_sparse_coding(U, D);
A = X*X';
B = U*X';
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, dict_iter);
norm(D*X - U, "fro")
num_iter
nnz(updated)
norm(D - D1, "fro")

X = lasso_sparse_coding(U, D);
A = X*X';
B = U*X';
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, dict_iter);
norm(D*X - U, "fro")
num_iter
nnz(updated)
norm(D - D1, "fro")
%{
X = lasso_sparse_coding(U, D);
A = X*X';
B = U*X';
[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, max_iter);
norm(D*X - U, "fro")
%}


%{
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
X = lasso_sparse_coding(U, D);
A = sum_of_outer_prods(X, X);
B = sum_of_outer_prods(U, X);
%}
%repeat update_dictionary()
%[D, D_prev, updated, num_iter] = dictionary_update(D, A, B, tol, max_iter);

%evaluate results of dictionary update  
%{
histogram(D)
max(max(D))
norm(D, "fro")
norm(D - D_prev, "fro")
nnz(updated)
num_iter


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




