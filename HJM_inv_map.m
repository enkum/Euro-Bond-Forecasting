% Copyright   : Ebenezer Nkum,Michael Pokojovy, and Thomas M. Fullerton, Jr. (2023)
% Version     : 1.0
% Last edited : 08/26/2023
% License     : Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
%               https://creativecommons.org/licenses/by-sa/4.0/

function [I_sigma_hat, lambda_hat, var_rat] = HJM_inv_map(T_grid, X_grid, Y_obs, r_obs, pc_var)
    
    dx = X_grid(2) - X_grid(1);
    dt = T_grid(2) - T_grid(1);

    % Predictor vectors
    pred = zeros(length(T_grid) - 1, length(X_grid));
    
    for i = 2:length(T_grid)
        f1  = Y_obs(i,     :);
        f0  = Y_obs(i - 1, :);
        %Af0 = [0 (f0(2:end) - f0(1:(end - 1)))/dx];
        Af1 = [0 (f1(2:end) - f1(1:(end - 1)))/dx];
        
        pred(i - 1, :) = (f1 - f0)/dt + Af1 - r_obs(i);
    end
    
    %% Estimate I sigma
    % Transform to D(A)
    for i = 1:(length(T_grid) - 1)
        pred(i, :) = A(pred(i, :));
    end
    
    S = cov(pred);
    [V, D] = eig(S);
    
    V = V/sqrt(dx); % normalize the discrete L^2 norm
    d = diag(D)*dx; % scale the eigenvalues
    
    var_rat = cumsum(flip(d))./sum(d);
    n_mode  = min(find(var_rat >= pc_var));
    
    sigma_hat  = diag(sqrt(d((end - n_mode + 1):end)*dt))*V(:, (end - n_mode + 1):end)'; % Note the normalization factor dt
    I_sigma_hat = zeros(size(sigma_hat));
    
    % Transform back to H
    for i = 1:n_mode
        I_sigma_hat(i, :) = I(sigma_hat(i, :));
    end   
    
    %% Estimate lambda
    loc = mean(pred, 1);
    
    % Transform to D(A)
    A_I_sigma_hat2 = zeros(size(I_sigma_hat));
    for i = 1:n_mode
        A_I_sigma_hat2(i, :) = A(I_sigma_hat(i, :).^2);
    end
    
    lambda_hat = 2*sqrtm(A_I_sigma_hat2*A_I_sigma_hat2')\(A_I_sigma_hat2*loc') - 1;
    
    function res = A(y)
        res = [0 (y(2:end) - y(1:(end - 1)))/dx];
    end

    function res = I(y)
        res = dx*cumsum(y);
    end
end
