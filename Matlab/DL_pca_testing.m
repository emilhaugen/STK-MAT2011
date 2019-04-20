%
% perform dictionary learning based on PCA
% map each data vector to its PCA representation and do DL on new data
% 

%% prepare
clear 

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";
SUBSET_LEN = 896; % should be > 64 to avoid too many zeros in diag(A)
BLOCK_LEN = 32; % data subset will have SUBSET_LEN^2 / BLOCK_LEN^2 vectors
%%   
% read data into block form
[U, restored] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); % < 1.5 sec


% set params for dictionary learning routine
CODE_LEN = 100; % Number of atoms 
lambda = .1; % training time increases as lambda decreases
dict_tol = 1e-6;
dict_iter = 20;
max_iter = 3;
%% PCA
% use eig.vectors of U*U' containing more than explain_tol of 
% total variance in data trace(U*U')
explain_tol = 1e-4;
[V, P, E, K] = PCA(U, explain_tol);

% map each data vector U(:,j) (each patch) to its low dim. PCA coefficients
U_pca = V' * U;
r = patches_to_original(V*U_pca, BLOCK_LEN, SUBSET_LEN, SUBSET_LEN);
 %% DictionaryLearning with PCA
[D, X, A, B, updated, n, i] = dictionary_learning(U_pca, CODE_LEN, ...
                        lambda, dict_tol, max_iter, dict_iter, restored, V);

 
%%
dict = V*D;
rec = patches_to_original(dict*X, BLOCK_LEN, SUBSET_LEN, SUBSET_LEN);
subplot(1, 2, 1)
imagesc(rec)
colorbar();
subplot(1, 2, 2)
imagesc(restored)
colorbar()
norm(restored - rec)
%% 
dict = V*D;
factors = make_dict_factors(dict, SUBSET_LEN, BLOCK_LEN); 
%% residuals covariance, and taper
orig_image = dlmread("Original.txt");
orig_image = orig_image(1:SUBSET_LEN, 1:SUBSET_LEN);
residual = orig_image - restored;
% estimate 
C = real(ifftn(abs(fftn(residual-mean(residual(:)))).^2))/numel(residual);
C(1,1) % these two 
var(residual(:)) % should be the same

L = fftn(C);

% make y*, D*
y_star = real(ifftn(sqrt(inv(L)) * fft(orig_image)));

imagesc(factors(:,:,1))
%% d_star
tic
for j=1:CODE_LEN % approx 1 min!
    factors(:,:,j) = real(ifftn(sqrt(inv(L)) * fft(factors(:,:,j))));
end
toc
imagesc(factors(:,:,1))
colorbar()
%%
size(D)







