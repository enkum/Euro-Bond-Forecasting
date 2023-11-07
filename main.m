% Copyright   : Ebenezer Nkum,Michael Pokojovy, and Thomas M. Fullerton, Jr. (2023)
% Version     : 1.0
% Last edited : 08/26/2023
% License     : Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
%               https://creativecommons.org/licenses/by-sa/4.0/

%% Set random seed for reproducibility
rng(1);

% Load and process 2018 data
yield2018 = readmatrix("yield_2018_euro_yield.csv");
 
yield2018 = yield2018(:, 3:(end - 1));

% Fill missing values using linear interpolation/ customize to your data
%search for where the data is empty/
% check the missing values in any of the columns from 3

T = 1 : 364;
dif = setdiff(T, find(isnan(yield2018(T , 3))));

% Get the position of the missing values
U = find(isnan(yield2018(T, 3)));

for p = 1:length(U)
    for k = 3:25
        yield2018(U(p), k)=interp1(dif, yield2018(dif,k), U(p),'linear','extrap');
    end 
end

% size of the data
[n, l] = size(yield2018);
month = [1 2 3 4 5 6 7 8 9 10 11 12*[1 2 3 4 5 6 7 8 9 10 15 20 25 30]];

I = 1:length(month);
I = setdiff(I, [1 2]);
for i = 1:n
    yield2018(i, 1) = interp1(month(I), yield2018(i, I), month(1),'linear','extrap');
    yield2018(i, 2) = interp1(month(I), yield2018(i, I), month(2),'linear','extrap');
end


% %% Loading 2019 data
yield2019 = readmatrix("yield_2019_euro_yield.csv");
yield2019 = yield2019(:, 3:(end - 1));


dif = setdiff(T,find(isnan(yield2019(T , 3))));
U = find(isnan(yield2019(T , 3)));
for p = 1 : length(U)
    for k = 3 : 25
        yield2019(U(p), k)=interp1(dif, yield2019(dif,k), U(p),'linear','extrap');
    end 
end
[q, h] = size(yield2019);
for i = 1:q
    yield2019(i, 1) = interp1(month(I), yield2019(i, I), month(1),'linear','extrap');
    yield2019(i, 2) = interp1(month(I), yield2019(i, I), month(2),'linear','extrap');
end


time_grid    = (1:n)*(12/n);
horizon_grid = month;

figure(1);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 28; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);

[Time, Horizon] = meshgrid(time_grid(:), horizon_grid);

surf(Time, Horizon, yield2018');
xlabel({'Calendar time $t$', '(in months)'}, 'FontSize', 25, 'interpreter', 'latex');
ylabel({'Time to maturity $x$', '(in months)'}, 'FontSize', 25, 'interpreter', 'latex');
zlabel({'Yield rate $y_{t}(x)$', '(in \%)'}, 'FontSize', 25, 'interpreter', 'latex');

n_month = 25;
month = month(1:n_month);

% Preparing the forward curves
T_grid     = linspace(0, 12, n);
X_grid     = linspace(0, month(end),  month(end)*361/ month(end));

dx = X_grid(2) - X_grid(1);
dt = T_grid(2) - T_grid(1);



yield_obs = zeros(length(T_grid), length(X_grid));
Y_obs     = zeros(size(yield_obs));

for i = 1:length(T_grid)
    yield_obs(i, :) = interp1([0 month], [yield2018(i, 1) yield2018(i, 1:n_month)], X_grid, 'spline');
    Y_obs(i, :)     = X_grid.*yield_obs(i, :);
end


yield_obs_19 = zeros(length(T_grid), length(X_grid));
Y_obs0     = zeros(size(yield_obs_19));

for i = 1:length(T_grid)
    yield_obs_19(i, :) = interp1([0 month], [yield2019(i, 1) yield2019(i, 1:n_month)], X_grid, 'spline');
    Y_obs0(i, :)     = X_grid.*yield_obs_19(i, :);
end



I_T = ceil(linspace(1, length(T_grid), 150));
I_X = ceil(linspace(1, length(X_grid), 100));


% Plotting integrated forward rates 

figure(2);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 28; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);

[Time, Horizon] = meshgrid(T_grid(I_T), X_grid(I_X));
surf(Time, Horizon, Y_obs(I_T, I_X)');

xlabel({'Calendar time $t$', '(in months)'}, 'FontSize', 25, 'interpreter', 'latex');
ylabel({'Time to maturity $x$', '(in months)'}, 'FontSize', 25, 'interpreter', 'latex');
zlabel({'Integrated forward rate $Y(t, x)$', '(in $\% \times \textrm{month}$)'}, 'FontSize', 25, 'interpreter', 'latex');

% Estimation of Vasicek's model
r_obs = yield_obs(:, 1);
[alpha_hat, r_ast_hat, beta_hat] = Vasicek_inv_map(T_grid, yield_obs(:, 1));
display(['Vasicek''s model: alpha-hat = ', num2str(alpha_hat), ...
         ', r*-hat = ', num2str(r_ast_hat), ', beta-hat = ', num2str(beta_hat)]);

% Inverse problem for the abstract Heath-Jarrow-Morton model
[I_sigma_HJM_hat, lambda_HJM_hat,var_rat] = HJM_inv_map(T_grid, X_grid, Y_obs, r_obs, 0.99);

display(['HJM model: lambda-hat = ', num2str(lambda_HJM_hat')]);

n_mode = size(I_sigma_HJM_hat, 1);

% Plot principal curves
figure(3);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 28; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);


hold on;

xlabel('$x$', 'FontSize', 25, 'interpreter', 'latex');
ylabel('Principal curves $\mathcal{I}_{x} \sigma_{n}$', 'FontSize', 25, 'interpreter', 'latex');

axis([min(X_grid) max(X_grid) 1.05*min(I_sigma_HJM_hat, [], 'all') 1.05*max(I_sigma_HJM_hat, [], 'all')]);

for i = 1:n_mode
    plot(X_grid, I_sigma_HJM_hat(i, :), 'LineWidth', 2, 'MarkerSize', 0.01);
end

legend({'$\mathcal{I}_{x} \sigma_{1}$', '$\mathcal{I}_{x} \sigma_{2}$', '$\mathcal{I}_{x} \sigma_{3}$','$\mathcal{I}_{x} \sigma_{4}$'}, ...
        'FontSize', 25, 'interpreter', 'latex', 'Location', 'NorthWest');

%% Prediction for Y
n_rep = 10000;
conf  = 0.998;

T_pred = 30; %Correspond to Jan 31, 2019
date   = 'January 31, 2019';


T_grid_pred = T_grid(1:T_pred);

t_pred = T_grid_pred(end);

% Predicting the short rates f_{t}(0) with Vasicek's model
pred0 = Vasicek_fwd_map(T_grid_pred, yield_obs_19(1,1), alpha_hat, r_ast_hat, beta_hat, n_rep);

% Predicting integrated forward rates with HJM model
pred  = HJM_fwd_map(T_grid_pred, X_grid, Y_obs0(1, :), pred0, I_sigma_HJM_hat, lambda_HJM_hat, n_rep);

yield_pred = zeros(length(X_grid), n_rep);

for j = 1:n_rep
    yield_pred(:, j) = pred(end, :, j)./X_grid;
    yield_pred(1, j) = yield_pred(2, j);
end

yield_lq  = zeros(size(X_grid));
yield_hq  = zeros(size(X_grid));
yield_avg = zeros(size(X_grid));

for i = 1:length(X_grid)
    yield_lq(i)  = quantile(yield_pred(i, :), (1 - conf)/2);
    yield_hq(i)  = quantile(yield_pred(i, :), 1 - (1 - conf)/2);
    yield_avg(i) = mean(yield_pred(i, :));
end


%%%%% calculating integrated mean square error %%%%%%

% Get the first 30 days for Jan 
% Define constants and parameters
n_rep = size(pred, 3);

num_time_steps = 30; % Number of time steps
num_grid_points = 361; % Number of grid points
% Assuming X_grid is a vector of grid points

% Extract relevant data from yield_obs_19_30
yield_obs_19_30 = yield_obs_19(1:num_time_steps, :);
obs_data = yield_obs_19_30(:, 2:num_grid_points);

% Initialize arrays for u and b
u = zeros(num_time_steps, 2);
b = zeros(num_time_steps, 1);

% Loop over time steps
for t = 1:num_time_steps
    % Calculate squared differences and sum them
    diff = (obs_data(t, :)) - (squeeze(pred(t, 2:num_grid_points, :))') ./ (X_grid(2:num_grid_points));
    b(t) = sum(diff(:).^2) * dx;

    % Normalize b
    b(t) = b(t)*1 / 10000;

    % Store results in u
    u(t, :) = [t, b(t)];
end


%% One-month prediction
figure(4);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 28; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);

xlabel({'Calendar Day in January 2019'}, 'FontSize', 15, 'interpreter', 'latex');
ylabel({'MISE'}, 'FontSize', 15, 'interpreter', 'latex');


axis([0 31 0 15]); % Modified y-axis limits

plot(u(:,1), u(:,2), 'k', 'LineWidth', 3);



% Find the indices 
X_grid_ext =X_grid;
I_X = find(ismember(X_grid_ext, X_grid));


%% One-month prediction
figure(5);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 36; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);



hold on;

xlabel({'Time to maturity $x$', '(in months)'}, 'FontSize', 15, 'interpreter', 'latex');
ylabel({['Predicted yields $y_{t}(x)$ on $t = \textrm{', date, '}$'], '(in \%)'}, 'FontSize', 15, 'interpreter', 'latex');


% axis([min(X_grid_ext(I_X)) max(X_grid_ext(I_X)) 0 4]);

axis([min(X_grid_ext(I_X)) max(X_grid_ext(I_X)) -1 4]); % Modified y-axis limits

%plot(1:25, yield2019(T_pred, 1:n_month), 'k', 'LineWidth', 3);

plot(X_grid, interp1(month, yield2019(T_pred, 1:n_month), X_grid,'spline'), 'k', 'LineWidth', 3);
% plot(X_grid, interp1(month, yield2019(T_pred1, 1:n_month), X_grid,'spline'), 'k', 'LineWidth', 3);
plot(X_grid_ext(I_X), yield_lq(I_X),  'k:',  'LineWidth', 2);
plot(X_grid_ext(I_X), yield_hq(I_X),  'k:',  'LineWidth', 2);
plot(X_grid_ext(I_X), yield_avg(I_X), 'k-.', 'LineWidth', 2);

for i = 1:5
    plot(X_grid_ext(I_X), yield_pred(I_X, i));
end

legend({'Observed yield curve',...
       ['Lower $', num2str(100*conf), '\%$ pointwise prediction bound'],...
       ['Upper $', num2str(100*conf), '\%$ pointwise prediction bound'],...
       'Estimated mean yield curve',...
       'Five sample yield curve forecasts'}, ...
        'FontSize', 10, 'interpreter', 'latex', 'Location', 'NorthWest');
  

%% One-month avg. rate curve estimation
figure(6);
set(gcf, 'PaperUnits', 'centimeters');
xSize = 36; ySize = 16;
xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [0 0 xSize*50 ySize*50]);


% Observed
subplot(1, 2, 1);
hold on;


xlabel({'Calendar time $t$', '(in months)'}, 'FontSize', 10, 'interpreter', 'latex');
ylabel({'Time to maturity $x$', '(in months)'}, 'FontSize', 10, 'interpreter', 'latex');
zlabel({'Obs. yield rate $y_{t}(x)$', '(in \%)'}, 'FontSize', 10, 'interpreter', 'latex');

n_days_jan = 30; % number of days in january 2019, Jan 1 is a holiday

[Time, Horizon] = meshgrid(T_grid_pred(1:n_days_jan), month);

surf(Time, Horizon, yield2019(1:n_days_jan, 1:n_month)');
axis([min(T_grid_pred(1:n_days_jan)) max(T_grid_pred(1:n_days_jan)), min(month) max(month), 0.0 2.5]);
view(-20, 50);



% Estimated mean
% Estimated mean
subplot(1, 2, 2);
hold on

xlabel({'Calendar time $t$', '(in months)'}, 'FontSize', 10, 'interpreter', 'latex');
ylabel({'Time to maturity $x$', '(in months)'}, 'FontSize', 10, 'interpreter', 'latex');
zlabel({'Est. mean yield rate $\widehat{\mathrm{E}\big[y_{t}(x)\big]}$', '(in \%)'}, 'FontSize', 10, 'interpreter', 'latex');


mean_yield_pred = mean(pred, 3);


for i = 1:size(pred, 1)
    mean_yield_pred(i, :) = mean_yield_pred(i, :)./X_grid_ext;
end



I_T = 1:n_days_jan;
I_X = ceil(linspace(1, length(X_grid_ext), 20));

[Time, Horizon] = meshgrid(T_grid_pred(I_T), X_grid_ext(I_X));

surf(Time, Horizon, mean_yield_pred(I_T, I_X)');
axis([min(T_grid_pred(I_T)) max(T_grid_pred(I_T)), min(X_grid_ext(I_X)) max(X_grid_ext(I_X)), 0.0 2.5]);
view(-20, 50);

% var_rat
%print('main_3','-depsc')