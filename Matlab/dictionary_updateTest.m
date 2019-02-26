%
%  Script to test dictionary_update function.
%  
%
%
%
%
%
SUBSET_LEN = 64;
BLOCK_LEN = 8;
CODE_LEN = 20; %length of sparse code vectors

U = ascii_to_data_matrix("Restoration.txt", SUBSET_LEN, BLOCK_LEN);
Udim = size(U);
m = Udim(1);
T = Udim(2); %using notation T for no. of data vectors as in paper

D = init_dict(U, CODE_LEN);
 
tol = 1e-6;
max_iter = 5;

X = zeros(CODE_LEN, T);
X = lasso_sparse_coding(U, D, X, T);
A = sum_of_outer_prods(X, X);
B = sum_of_outer_prods(U, X);

size(D)






