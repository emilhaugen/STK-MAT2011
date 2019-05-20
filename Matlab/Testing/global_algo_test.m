%
%  Global script to test convergence.
%
%set constants

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 256; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 8; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; % length of sparse code vectors, no. of colums in dictionary

% convergence parameters
dict_update_tol = 1e-4;
dict_iter = 50; % max iterations in dict. update routine
glob_iter = 7; % max iterations in global algorithm

% read ascii data
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN);
T = length(U(1,:));
% randomly initialize dictionary
%D = init_dict(U, CODE_LEN, false);
%{
tic
glob_updated_col = zeros(CODE_LEN, 1);
for j = 1:glob_iter
    fprintf("\nGLOB ITER %d\n", j)
    [X, A, B] = lasso_sparse_coding(U, D);
    diff = norm(U - D*X, "fro");
    fprintf("U - DX=%0.5e\n", diff);
    i = 0;
    [D, D_prev, updated] = dictionary_update_util(D, A, B);
    while i < dict_iter && norm(D-D_prev, "fro") > dict_update_tol
        [D, D_prev, updated] = dictionary_update_util(D, A, B);
        glob_updated_col = glob_updated_col + updated;
        i = i + 1;
    end
    fprintf("%d dictionary update iters\n", i)
end
[X, A, B] = lasso_sparse_coding(U, D);
diff = norm(U - D*X, "fro");
fprintf("U - DX=%0.5e\n", diff);
toc
%}












