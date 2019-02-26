%
%  Script to test dictionary_update function.
%  
%
%
%
%
%
addpath("Utility_Functions");
addpath("Preprocess");


SUBSET_LEN = 64;
BLOCK_LEN = 8; %data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; %length of sparse code vectors, no. of colums in dictionary

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

if abs(min(diag(A))) < 1e-9
    error("Error: Zero diagonal element in A")
end 
%%{
for i = 1:1
    D_prev = D;
    for j = 1:CODE_LEN %initial update
        dj = B(:,j) - D*A(:,j) + D(:,j)*A(j,j); 
        D(:,j) = dj / (norm(dj) * A(j,j));
    end 
end    

   
%{
norm(D - D_prev, "fro")
max(max(D))
max(max(D_prev))
%}
z = zeros(T);
z = z(:, 1);
for p = 1:T
    z(p) = CODE_LEN - nnz(X(:,p));
end    
z    
%%}




