addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 8; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
CODE_LEN = 100; % length of sparse code vectors, no. of colums in dict.
   
% read data into block form
[U, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % time < 1.5 sec
T = length(U(1,:));

% Principal Component Analysis (slow by-hand method, use SVD for speed)
M = mean(U, 2);
V = U - M; % put data on mean deviation form

[P, E]=eig(V*V');

% will use principal components whose corresponding ratio 
% eigenvalue/trace(E) is at least explian_tol
explain_tol = 1e-3; 

% get no. of princ. comp. that have eigenvalue above threshold
pc_filter = diag(E) / trace(E) > explain_tol;
num_pc = sum(pc_filter); % will be length of sparse codes

fprintf("%d significant PC's of %d\n", num_pc, length(pc_filter))

% use principal components corresponding to num_pc highest eigenvalues
D = P(:,end-num_pc+1:end);

pca_coefficients = D' * V;
reconstruct = patches_to_original(D * pca_coefficients);

diff = norm(original - reconstruct, "fro");

fig = imagesc(original);
colorbar
title("Original data")
saveas(fig, "Original.png")

fig = imagesc(reconstruct);
colorbar
title(sprintf("Reconstructed from %d PC's at tolerance %0.1e with %dx%d patches\nNorm(original-reconstruct, 'fro')= %0.2e", ...
                num_pc, explain_tol, BLOCK_LEN, BLOCK_LEN, diff))
saveas(fig, "PCA_Reconstruction.png")            




