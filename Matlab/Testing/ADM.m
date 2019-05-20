function [beta_null, i] = ADM(y, D, lambda, rho, rho_scale, tol, max_iter)
    % ADM implement Alternating Direction Method to find sparse
    %   coefficients in prediction by relaxation.
    %
    %   param y (vector): target       
    %
    %   param D (matrix): constant
    %
    %   param beta (vector): initial coefficients, usually just zero
    %   
    %   param lambda (float): fixed regularization parameter
    %
    %   param rho (float): relaxation parameter. Initially small, 
    %                       increases as algorithm iterates.
    %
    %   param rho_scale (float): mulitply rho by this each iteration 
    %
    %   param tol (float): error tolerance
    %   
    %   param max_iter (int): max. no. of ADM iterations
    %
    %   return beta (float): sparse coefficients 
    
    % initialize relaxation parameter to zero
    num_atoms = length(D(1,:));
    beta_null = zeros(num_atoms, 1); 
    DTD = D'*D;
    % initial estimate for beta
    beta = inv(DTD/rho + eye(num_atoms)) * (D'*y/rho + beta_null);
    for j = 1:num_atoms
        beta_null(j) = sign(beta(j)) * max(0, abs(beta(j))-lambda/rho);
    end
    i = 0;
    error = norm(y-D*beta)/norm(y);
    while error > tol && i < max_iter % repeat
      % take beta_null constant  
      %norm_beta = norm(beta)
      beta = inv(DTD/rho + eye(num_atoms)) * (D'*y/rho + beta_null);
      % taking beta constant, optimize wrt. beta_null:
      % using component wise rule,
      % beta_null = sgn(beta)*max(0, abs(beta) - lambda/rho)
      for j = 1:num_atoms
        beta_null(j) = sign(beta(j)) * max(0, abs(beta(j))-lambda/rho);
      end
      % take beta_null constant
      %error = norm(y-D*beta)/norm(y);
      rho = rho_scale * rho;
      i = i + 1;
    end
end

