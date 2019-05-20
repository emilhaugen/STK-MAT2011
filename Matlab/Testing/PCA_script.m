addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 32; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
   
% read data into block form
[U, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % time < 1.5 sec
T = length(U(1,:));

% will use principal components whose corresponding ratio 
% eigenvalue/trace(E) is at least explian_tol
explain_tol = 1e-4; % epsilon

[V, P, E, U] = PCA(U, explain_tol);

y = diag(E);
k = length(y(y > 1e+3))
k = 200;
y_big = y(end-k+1:end);
n = length(y_big);
fig = scatter(linspace(1, n, n), y_big);
title(sprintf("Biggest %d eigenvalues", n));
set(gca, "yscale", "log")
saveas(fig, "Plots/PcaPlots/Biggest.png")

y_small = y(1:end-k);
n = length(y_small);
fig = scatter(linspace(1, n, n), y_small);
title(sprintf("Smallest %d eigenvalues", n));
set(gca, "yscale", "log")
saveas(fig, "Plots/PcaPlots/Smallest.png")

%{
fig = imagesc(original);
colorbar
title("Original data (centered)")
saveas(fig, "Plots/Original.png")

fig = imagesc(reconstruct);
colorbar
title(sprintf("Reconstructed from %d PC's at tolerance %0.1e with %dx%d patches\nNorm(original-reconstruct, 'fro')= %0.2e", ...
                num_pc, explain_tol, BLOCK_LEN, BLOCK_LEN, diff))
            
saveas(fig, sprintf("Plots/PCA_Reconstruction-%dPCs-BLOCK_LEN=%d.png", num_pc, BLOCK_LEN))            
%}



