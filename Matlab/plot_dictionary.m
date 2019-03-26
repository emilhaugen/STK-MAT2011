function plot_dictionary(D, V, BLOCK_LEN, directory)
    % PLOT_DICTIONARY plot atoms of dictionary and reconstruction
    % 
    %   param D (matrix): dictionary
    %
    %   param V (matrix): eigenvectors to map PCA dictionary back to 
    %                      original format. If already on original format,
    %                      use identity matrix.
    %   
    %   param BLOCK_LEN (int): side length of square patches needed for
    %                           reconstruction
    %
    %   param directory (string): location to save plots
    % 
    %   return NOTHING
    
    if ~isfolder(directory)
        error(strcat("Could not find directory named ", directory));
    end    
    
    delete(strcat(directory, "*.png"));
    
    num_atoms = length(D(1,:)); % No. of atoms/columns in dictionary
    rec_dict = V*D; % Reconstructed dictionary
       
    for j = 1:num_atoms
        patch = reshape(rec_dict(:,j), BLOCK_LEN, BLOCK_LEN)';
        fig = imagesc(patch);
        colorbar;
        title(sprintf("Atom no. %d of %d\n", j, num_atoms));
        saveas(fig, strcat(directory, sprintf("Atom%03d.png", j)));
    end    
end

