addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 32; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; % length of sparse code vectors, no. of colums in dict.
   
% read data into block form
U = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % time < 1.5 sec
T = length(U(1,:));
% Principal Component Analysis (slow by-hand method, use SVD for speed)
M = mean(U, 2);
V = U - M; % put data on mean deviation form
[P, E]=eig(V*V');

% will use principal components whose corresponding ratio 
% eig.value/trace(D) is at least explian_tol
explain_tol = 1e-4; 

% get no. of princ. comp. that have eigenvalue above threshold
pc_filter = diag(E) / trace(E) > explain_tol;
num_pc = sum(pc_filter) % will be length of sparse codes

% use principal components corresponding to num_pc highest eigenvalues
D = P(:,end-num_pc+1:end); % dictionary matrix

% use LASSO to find sparse code coefficients X^1, ..., X^T
tic
[X, A, B] = lasso_sparse_coding(V, D);
toc
%3000 fro


