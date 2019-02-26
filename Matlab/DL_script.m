%
%  High level script to test and debug functions.    
%
SUBSET_LEN = 128;
BLOCK_LEN = 8;
CODE_LEN = 50; %length of sparse code vectors

U = ascii_to_data_matrix("Restoration.txt", SUBSET_LEN, BLOCK_LEN);


Udim = size(U);
m = Udim(1);
T = Udim(2); %using notation T for no. of data vectors as in paper

D = init_dict(U, CODE_LEN);

[Dict, X, A, B] = dictionary_learning(U, D, 1e-6, 20);

A
X




















