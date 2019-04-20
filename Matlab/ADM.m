function beta = ADM(y, D, lambda, rho, rho_scale, tol, max_iter)
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
    %                       increase exponentially as algorithm iterates.
    %
    %   param rho_scale (float): mulitply rho by this each iteration 
    %
    %   param tol (float): error tolerance
    %   
    %   param max_iter (int): max. no. of ADM iterations
    %
    %   return beta (float): sparse coefficients 
    
    % initialize relaxation parameter to zero
    n = length(D(:,1));
    beta_null = zeros(n); 
    
    % initial estimate for beta
    beta = inv(D' * D/rho + eye(n)) * (D'*y/rho + beta_null);
    i = 0;
    while norm(y-D*beta) + lambda(norm(beta_null, 1)) + ...
      rho*norm(beta - beta_null) > tol && i < max_iter
      % taking beta constant, optimize wrt. beta_null:
      % using component wise rule,
      % beta_null = sgn(beta)*max(0, abs(beta) - lambda/rho)
      for j = 1:n
        beta_null(j) = sign(beta(j)) * max(0, abs(beta(j))-lambda/rho);
      end
      beta = inv(D'*D/rho + eye(n)) * (D'*y/rho + beta_null);
      rho = rho*rho_scale;
      i = i + 1;
    end
end

