% Copyright   : Ebenezer Nkum,Michael Pokojovy, and Thomas M. Fullerton, Jr. (2023)
% Last edited : 08/26/2023
% License     : Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
%               https://creativecommons.org/licenses/by-sa/4.0/

function sol = Vasicek_fwd_map(T_grid, y0, alpha, r_ast, beta, n_rep)
    dt = T_grid(2) - T_grid(1);

    sol = zeros(length(T_grid), n_rep);
    
    a = alpha*dt;
    b = alpha*r_ast*dt;
    c = beta*sqrt(dt);
    
    for k = 1:n_rep 
        sol(1, k) = y0;
        for j = 2:length(T_grid)

            sol(j, k) = b + (1 - a)*sol(j - 1, k) + c*randn(1);
           
        end
    end
end
