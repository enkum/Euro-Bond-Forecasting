% Copyright   : Ebenezer Nkum,Michael Pokojovy, and Thomas M. Fullerton, Jr. (2023)
% Version     : 1.0
% Last edited : 08/26/2023
% License     : Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
%               https://creativecommons.org/licenses/by-sa/4.0/

function [alpha_hat, r_ast_hat, beta_hat] = Vasicek_inv_map(T_grid, r_obs)
    dt = T_grid(2) - T_grid(1);
    
    r_cur = r_obs(2:end);
    r_lag = r_obs(1:(end - 1));
    
    m = length(T_grid) - 1;
    
    alpha_hat = (-1/dt)*log((m*sum(r_cur.*r_lag) - sum(r_cur)*sum(r_lag))/(m*sum(r_lag.^2) - (sum(r_lag))^2));
    
    expa  = exp(alpha_hat*dt);
    expai = 1.0/expa;
    
    r_ast_hat = 1/(m*(1 - expai))*(sum(r_cur) - expai*sum(r_lag));
    
    beta_hat = sqrt(2*alpha_hat/(m*(1 - expai^2))*sum((r_cur - r_lag*expai - r_ast_hat*(1 - expai)).^2));
end
