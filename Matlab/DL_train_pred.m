%% prepare
clear 

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");


%% set constants and load data

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 32; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors

[U, restored] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % < 1.5 sec

% set params for dictionary learning routine
CODE_LEN = 175; % Number of atoms 
lambda = .05; % training time increases as lambda decreases
dict_tol = 1e-6;
dict_iter = 40;
max_iter = 20;

%% PCA
% use eig.vectors of U*U' containing more than explain_tol of 
% total variance in data trace(U*U')
explain_tol = .51e-4;
[V, P, E, K] = PCA(U, explain_tol);

% map each data vector U(:,j) (each patch) to its low dim. PCA coefficients
U_pca = V' * U;

 %% DL on PCA data
tic
[D, X, A, B, errors, num_iter] = dictionary_learning(U_pca, CODE_LEN, ...
                        lambda, dict_tol, max_iter, dict_iter, ...
                        restored, V, SUBSET_LEN, BLOCK_LEN);
stop_time = toc        
% rescale dictionary
dictionary = V*D;

%% residuals and compute sigma
original = dlmread("Original.txt");
original = original(1:SUBSET_LEN, 1:SUBSET_LEN);

residual = original - restored;
% estimate covariance of residuals
C = real(ifftn(abs(fftn(residual-mean(residual(:)))).^2))/numel(residual);

s = sigma_matrix(C', BLOCK_LEN);
imagesc(s);

L = (s^(-0.5));

%% predict to recover signal from residual

lambda = 0.015;
recovered = lasso_predict(residual, dictionary, L, lambda, BLOCK_LEN);


%% display recovered signal and updated estimate along with resotored data

imagesc(restored); colorbar
caxis([-5, 5])
title("Initial estimate");
figure()

imagesc(recovered); colorbar
caxis([-5, 5])
title("Recovered signal") 
figure()
imagesc(recovered + restored)
title("Updated estimate")
colorbar()
caxis([-5, 5])


%% zoom in on predictions

% add recovered signal to initial estimate 'restored'
u_est = recovered + restored; 

rows = 700:896;
cols = 370:550;

imagesc(restored(rows, cols));
title("initial estimate")
caxis([-6, 6])
colorbar
figure()
imagesc(u_est(rows, cols));
title("Restored + recovered signal")
caxis([-6, 6]);
colorbar






