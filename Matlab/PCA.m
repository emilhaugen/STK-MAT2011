function [signf_comp, P, D, U]  = PCA(U, tol)
    % PCA Perform simple principle component analysis
    %   
    %   param U (matrix): Input data matrix
    % 
    %   param tol (float): Tolerance for significance of eigenvectors.
    %     Eigenvectors of U*U' whose corresponding eigenvalue divided 
    %     by the sum of all eigenvalues exceeds tol will be used as 
    %     signinficant principal components of the data.
    %                       
    %   return signf_comp (matrix): significant components
    %
    %   return P (matrix): all components    
    %
    %   return D (matrix): diagonal matrix with ALL eigenvalues
    %
    %   return U (matrix): centered data
    %
    
    if norm(mean(U, 2)) > 1e-14 % put on mean deviation form if not already
        fprintf("Centering data\n");
        U = U - mean(U, 2);
    end 
    [P, D] = eig(U*U');
    
    % number of significant components
    num_signf_pc = sum(diag(D) / trace(D) > tol);
    
    if num_signf_pc > 0   
        fprintf("%d significant PCs\n", num_signf_pc)
        signf_comp = P(:,end-num_signf_pc+1:end); 
    else     
        error(sprintf("No significant principal components at tolreance %.1e", tol));  
end

