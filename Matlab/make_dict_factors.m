function [D_factors] = make_dict_factors(D, image_length, patch_len)
    % MAKE_DICT_FACTORS create dictionary factors used in filtering
    %   One factor for each atom. Factors are square matrices 
    %    of same dimension as square image. 
    num_atoms = length(D(1,:));
    D_factors = zeros(image_length, image_length, num_atoms);
    n = patch_len / 2; % PATCH_LEN = BLOCK_LEN

    for k = 1:num_atoms
      patch = reshape(D(:,k), patch_len, patch_len)';
      % subdivide patch
      top_right = patch(1 : n, n + 1 : end);
      bottom_right = patch(n + 1 : end, n + 1 : end);
      bottom_left = patch(n + 1 : end, 1 : n);
      top_left = patch(1 : n, 1 : n);
      % map patch subsets to factor corners
      N = image_length - n + 1;
      D_factors(N : end, 1 : n, k) = top_right; % bottom left
      D_factors(1 : n, 1 : n, k) = bottom_right; % top left
      D_factors(1 : n, N : end, k) = bottom_left; % top right
      D_factors(N : end, N: end, k) = top_left; % bottom right
    end    
end

