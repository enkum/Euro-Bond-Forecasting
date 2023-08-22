% Copyright   : Michael Pokojovy & Valerii Maltsev (2020)
% Version     : 2.0
% Last edited : 11/27/2020
% License     : Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
%               https://creativecommons.org/licenses/by-sa/4.0/

function sol = HJM_fwd_map(T_grid, X_grid, Y0, y0, I_sigma, lambda, n_rep)
    k = size(I_sigma, 1);
    
    dt = T_grid(2) - T_grid(1);
    dx = X_grid(2) - X_grid(1);
    
    n = length(X_grid) - 1;
    
    A   = speye(n)/dx - sparse(2:n, 1:(n - 1), ones(1, n - 1), n, n)/dx;
    Adt = speye(n) + A*dt;
    
    sol = zeros(length(T_grid), length(X_grid), n_rep);
    
    for j = 1:n_rep
        sol(1, :, j) = Y0;
        
         for i = 2:length(T_grid)
%         dt = d(i)*(T_grid(2) - T_grid(1));
%             if mod(i,5)==1
%                 dt = T_grid(4) - T_grid(1);
%             else
%                 dt = T_grid(2) - T_grid(1); 
%             end
            dW = randn(1, k);
            
            % Solve with (semi-)implicit Euler-Maruyama
            rhs = sol(i - 1, 2:end, j) + ...
                  dt*(0.5*(1 + lambda')*(I_sigma(:, 2:end).^2) + y0(i - 1, j)) + sqrt(dt)*dW*I_sigma(:, 2:end);
            rhs = rhs';
            
            sol(i, :, j) = [0 (Adt\rhs)'];
        end
    end
end