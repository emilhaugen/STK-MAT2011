%
% perform dictionary learning based on PCA
% 

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 32; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
   
% read data into block form
[U, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % time < 1.5 sec

[V, P, E, U] = PCA(U, 1e-3); % explain_tol can be arbitrarily low here

% use MAX_CODE_LEN or all components if fewer than MAX_CODE_LEN
% take MAX_CODE_LEN to be twice the number of significant components
MAX_CODE_LEN = 2 * length(V(1,:))
CODE_LEN = min(MAX_CODE_LEN, length(P(1,:))-1);
D = P(:, end - CODE_LEN + 1 : end); % initial dictionary

[X, A, B] = lasso_sparse_coding(U, D);
norm(U - D*X, "fro")

[D, D_prev, upd, i] = dictionary_update(D, A, B, 1e-8, 50);

norm(U - D*X, "fro")

[X, A, B] = lasso_sparse_coding(U, D);
norm(U - D*X, "fro")

[D, D_prev, upd, i] = dictionary_update(D, A, B, 1e-8, 50);

norm(U - D*X, "fro")




 