
% Script to test that reading data from file with 
%   ascii_to_data_matrix() then using patches_to_original() on patched
%   output from first function yields original data, i.e. functions are 
%   inverses.

addpath("Help_Functions"); %make help functions available
addpath("Preprocess");
addpath("Data");

FILENAME = "Restoration.txt";

SUBSET_LEN = 896;
BLOCK_LEN = 32; % must divide SUBSET_LEN
   
% read data into block form
[patches, original] = ascii_to_data_matrix(FILENAME, SUBSET_LEN, BLOCK_LEN); 

% restore original format
original_from_patches = patches_to_original(patches);

norm(original_from_patches - original, "fro")



