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
 %% DictionaryLearning on PCA data
tic
[D, X, A, B, errors, num_iter] = dictionary_learning(U_pca, CODE_LEN, ...
                        lambda, dict_tol, max_iter, dict_iter, ...
                        restored, V, SUBSET_LEN, BLOCK_LEN);
stop_time = toc        
% rescale
dict = V*D;
%% errorplot
error_save_name = sprintf("Plots/Errorplots/numAtoms=%d-numPC=%iter=%d.png", CODE_LEN, length(V(1,:)), num_iter);
if isfile(error_save_name)
    delete(error_save_name);
end

err = errors(errors > 1e-10);
n = length(err);
k = 2;
t = linspace(k, n, n-k+1)
fig = plot(t, err(k:end));
xticks(t(mod(t, 2)==0));
s1 = sprintf("Relative error = %.5f after %d iterations with lambda=%.3f.", err(end), n, lambda);
s2 = sprintf("\nUsing %d atoms and %d significant principal components.", CODE_LEN, length(V(1,:)));
s3 = sprintf("\nTotal training time: %.1f seconds", stop_time);
title(strcat(s1, s2, s3));
xlabel("Number of iterations");
ylabel("Relative error");
saveas(fig, error_save_name);

%% reconstructplot
rec_save_name = sprintf("Plots/ReconstructPlots/rec-numAtoms=%d-numPC=%d.png", CODE_LEN, length(V(1,:)));
if isfile(rec_save_name)
    delete(rec_save_name);
end

dict = V*D;
rec = patches_to_original(dict*X, BLOCK_LEN, SUBSET_LEN, SUBSET_LEN);
imagesc(restored)
title("Estimate")
colorbar
t = 16;
k = 16;
q = 28;
start = t*32+1:t*32+32;
stop = k*32+1:k*32+32;
imagesc(rec(start, stop));
title("Estimate");
figure()
imagesc(restored(start, stop))
title("Training data")
figure()
a = t*q + k + 1
imagesc(reshape(U(:,a), 32, 32)')
%% residuals covariance, and taper
orig_image = dlmread("Original.txt");
orig_image = orig_image(1:SUBSET_LEN, 1:SUBSET_LEN);
rest = dlmread("Restoration.txt");
rest = rest(1:SUBSET_LEN, 1:SUBSET_LEN);
residual = orig_image - rest;
% estimate covariance of residuals
C = real(ifftn(abs(fftn(residual-mean(residual(:)))).^2))/numel(residual);
C(1,1) % these two 
var(residual(:)) % should be the same 

B = (residual-mean(residual(:)).^2)/numel(residual);

imagesc(B);

%% make sigma matrix for residual patches in ADM
s = sigma_matrix(C, BLOCK_LEN);


imagesc(s);


%% tapering, lambda^-1/2 

xV = cat(2, [1:448], flip([1:448]));
gT=exp(-0.5*(xV'/20).^2)*exp(-0.5*(xV/200).^2); % gaussian taper
imagesc(fftshift(gT))
C_taper = C.*gT;
epsilon = 0.005622;
L_half_inv = 1;%./sqrt(abs(fftn(C_taper))+epsilon);
imagesc(fftshift(L_half_inv))

%% d_star-factors, y* and D matrix
y_star = real(ifftn(L_half_inv.*fftn(orig_image))); % target
imagesc(y_star);

%% d_star-factors, and D matrix
image_len = numel(y_star);
num_atoms = CODE_LEN;
D = zeros(image_len, num_atoms);
factors = make_dict_factors(dict, SUBSET_LEN, BLOCK_LEN); 
tic
for j=1:num_atoms % 14 sec
    tmp = real(ifftn(fftn(factors(:,:,j)).*L_half_inv)); % d_i*
    factors(:,:,j) = tmp;
    D(:,j) = tmp(:);
end
toc
%% inspect factors
imagesc(fftshift(factors(:,:,1)))
colorbar();
for i=1:num_atoms
  tmp = ifftn(fftn( factors(:,:,i)).*fftn(y_star));   
  dsy(:,:,i)= tmp;
  D_new(:,i) = tmp(:);
end
imagesc(dsy(:,:,1))


%% estimate y* with D*X by ADM algo
[beta, niter] = ADM(y_star(:), D, 1e-6, 1.1, 1.5, 1e-4, 20);

%% single patch
k = 9;
t = 13;
row_num = t*32+1:t*32+32;
col_num = k*32+1:k*32+32;
patch = residual(row_num, col_num);
beta = ADM(L*patch(:), L*dict, 0.1, 0.1, 1.2, 1e-4, 100);
beta_l = lasso(L*dict, L*patch(:), "Lambda", 0.0625)

imagesc(reshape(dict*beta_l, 32, 32))
title("lasso")
figure()
imagesc(reshape(dict*beta, 32, 32))
figure()
imagesc(patch)
title("patch")

%% patch wise
reconstruct = zeros(896, 896);
s = sigma_half_inv;
target = residual;
rest_patch = rest(800:831, 390:(390+31));
patch = residual(800:831, 390:(390+31));
imagesc(rest_patch)
figure()
imagesc(patch)

beta = ADM((s^(-0.5))*patch(:), (s^(-0.5))*dict, 1, 0.1, 1.2, 1e-4, 30)
imagesc(reshape(dict*beta, 32, 32)); colorbar
%%
s = sigma_matrix(C, BLOCK_LEN);
L=(s^(-0.5));
target = residual;
%%
%reconstruct = zeros(896, 896);
lasso_reconstruct = zeros(896, 896);
lambda = .01;
rho = 0.1;
rho_scale = 1.2;
tic
for t = 0:27
    t
    for k = 0:27
        row_num = t*32+1:t*32+32;
        col_num = k*32+1:k*32+32;
        patch = target(row_num, col_num);
        %beta = ADM(L*patch(:), L*dict, lambda, 0.01, 1.2, 1e-4, 100);
        beta_lasso = lasso(L*dict, L*patch(:), "Lambda", 0.015);
        %rec = reshape(dict*beta, 32, 32);
        %reconstruct(row_num, col_num) = rec;
        lasso_reconstruct(row_num, col_num) = reshape(dict*beta_lasso, 32, 32);
    end
end  
toc
%%
b = 8;
start = start;
stop = stop;
imagesc(target(start, stop))
figure()
imagesc(reshape(reconstruct(start, stop), 32, 32))
%% inspect
imagesc(reconstruct); colorbar
caxis([-5, 5])
title(sprintf("rec:lambda=%.2f\nrho=%.2f-rscale=%.2f", lambda, rho, rho_scale))
figure()
imagesc(target); colorbar
caxis([-5, 5])
title("original") 
err = norm(target - reconstruct, "fro") / norm(target, "fro");
figure()
imagesc(target - reconstruct)
title(sprintf("diff:relative difference: %.5f", err));
figure()
imagesc(lasso_reconstruct)
title("Recovered signal")
colorbar()
caxis([-5, 5])

%%
a = 700:896;
b = 370:550;
imagesc(restored(a, b));
title("Restored")
caxis([-6, 6]);
colorbar
figure()
u_est = lasso_reconstruct + restored;
imagesc(u_est(a, b));
title("Restored + recovered signal")
caxis([-6, 6]);
colorbar
figure()
imagesc(orig_image)
title("orig")
caxis([-6, 6])
%%
% make y*, D*
y_star = real(ifftn(L_half_inv.*fftn(orig_image)));
y_starO = real(ifftn(L_half_inv.*fftn(restored)));
imagesc(y_star-y_starO); caxis([-1 1])

% make y*, D*
imagesc(y_star)
colorbar
imagesc(orig_image)
%%
xV = [0:447-448:-1]; %xv mxn then gt is nxn
gT=exp(-0.5*(xV'/20).^2)*exp(-0.5*(xV/200).^2)
imagesc(fftshift(gT))
imagesc(fftshift(C))
imagesc(fftshift(C.*gT))
imagesc(factors(:,:,1))
imagesc(fftshift(abs(fftn(C))))

imagesc(fftshift(abs(fftn(C.*gT))))

imagesc(residual);
title("Residual data");
colorbar

%% d_star-vectors
factors = make_dict_factors(dict, SUBSET_LEN, BLOCK_LEN); 
tic
for j=1:CODE_LEN % 14 sec
    factors(:,:,j) = real(ifftn(fftn(factors(:,:,j)).*L_half_inv ));
end
toc
%%
imagesc(factors(:,:,1))
colorbar()
%% analyze coeff.

p = length(X(1,:));

w = sum(abs(X) > 1e-5, 1); % w(i)=no. of atoms used by u(i) est. i < T
median(w)
fig = histogram(w, 35)
%s1 = sprintf("Histogram of number of atoms used in reconstruction of \n each of the %d patches.", p);
%s2 = sprintf("\nMedian: %d out of %d atoms", median(w), num_atoms);
%title(strcat(s1, s2));
xlabel("Number of atoms used in reconstruction of blocks");
ylabel("Frequency");
saveas(fig, "Plots/Coefficients/atoms-per-patch.png");
figure()
w = sum(abs(X) > 1e-5, 2); % w(i)=no. of times atom (i) is used, i < k
fig = histogram(w, 30)
%s1 = sprintf("Histogram of how many reconstructions each atom occurs in \n out of %d reconstructed patches", p);
%s2 = sprintf("\nMedian: %d out of %d reconstructions", median(w), p);
%title(strcat(s1, s2));
xlabel("Number of times each atom is used in reconstruction");
ylabel("Frequency");
saveas(fig, "Plots/Coefficients/patches-per-atom.png");





