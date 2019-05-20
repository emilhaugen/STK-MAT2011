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
[U, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); 

data = patches_to_original(U, BLOCK_LEN, SUBSET_LEN, SUBSET_LEN);

[U_orig, original] = ascii_to_data_matrix("Original.txt", SUBSET_LEN, BLOCK_LEN);
imagesc(original)
imagesc(data)
residual = original - data;
imagesc(residual)
data = residual;
C=real(ifftn(abs(fftn(data-mean(data(:)))).^2))/numel(data);
C(1,1) % these two 
var(data(:)) % should be same

imagesc(C)
imagesc(fftshift(C)); colorbar

subplot(1,2,1)
plot(ifftshift(C(1,:)))
subplot(1,2,2)
plot(ifftshift(C(:,1)))





