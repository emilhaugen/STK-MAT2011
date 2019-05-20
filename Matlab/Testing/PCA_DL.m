%
% perform dictionary learning based on PCA
% map each data vector to its PCA representation and do DL on new data
% 
clear 

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 16; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
   
% read data into block form
[U, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % < 1.5 sec

% use eig.vectors of U*U' containing more than explain_tol of 
% total variance in data trace(U*U')
explain_tol = 1e-3;
[V, P, E, K] = PCA(U, explain_tol);

% map each data vector U(:,j) to its low dim. PCA coefficients
U_pca = V' * U;

% set params for dictionary learning routine
CODE_LEN = 150; % Number of atoms 
lambda = 0.05; % training time increases as lambda decreases
dict_tol = 1e-2;
dict_iter = 20;
max_iter = 20;

[D, X, A, B, updated, n, i] = dictionary_learning(U_pca, CODE_LEN, ...
                                lambda, dict_tol, max_iter, dict_iter);

write_dir = strcat("Plots/PcaDL/", ...
                    sprintf("numAtoms=%d-lambda=%.2e/", CODE_LEN, lambda));                            
  
if ~isfolder(write_dir)
    mkdir(write_dir);               
end    
                
plot_dictionary(D, eye(length(D(:,1))), BLOCK_LEN, write_dir);

diff = norm(U - D*X, "fro");
rel_diff = diff / norm(U, "fro");
fig = imagesc(patches_to_original(D*X, BLOCK_LEN, SUBSET_LEN, SUBSET_LEN));
title(strcat(sprintf("norm(original - reconstruct) = %.3e\n", diff), ...
                sprintf("relative norm = %.3e\n", rel_diff), ...
                sprintf("After %d DL iterations\n", i)));
colorbar;            
saveas(fig, strcat(write_dir, "reconstruct.png"));
                                        


%{
[X, A, B] = lasso_sparse_coding(U_pca, D, lambda);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);


[D, D_prev, updated, i] = dictionary_update(D, A, B, 1e-6, 10);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);

[X, A, B] = lasso_sparse_coding(U_pca, D, lambda);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);

[D, D_prev, updated, i] = dictionary_update(D, A, B, 1e-6, 10);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);

[X, A, B] = lasso_sparse_coding(U_pca, D, lambda);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);

[D, D_prev, updated, i] = dictionary_update(D, A, B, 1e-6, 10);
nrm = norm(U_pca, "fro");
diff = norm(U_pca - D*X, "fro"); 
fprintf("Abs.diff = %.3e, Rel.diff = %.3e\n", diff, diff/nrm);

manually try some lasso for a few of the data vectors
j = 230;
dataj = U_pca(:,j);
tic
[B, fitinfo] = lasso(D, dataj, "CV", 5);
toc
lassoPlot(B, fitinfo, "PlotType", "CV");
legend("show")
figure()
[mse, i] = min(fitinfo.MSE);
coef = B(:,i);
D * coef
U_pca(:,j);
r = D*coef;
i
mse
plot(dataj)
hold on 
plot(r)
hold on 
plot(dataj - r)
legend("orig", "rec", "diff")
hold off
%}
            