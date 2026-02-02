% Define your 2D data points (x, y, z)
% Replace these arrays with your own data
clear;
close all;
load('datasets/dataset2.mat');

F = @(x,y) sin(x).*cos(y);

x=vertices_0(1,:);
y=vertices_0(2,:);
z=vertices_0(3,:);

% Define the order of the B-spline in x and y directions
px = 2; % Order in x-direction (quadratic)
py = 2; % Order in y-direction (quadratic)

% Define the knot vectors for x and y
knots_x = [-3,-3 -3, -1,1, 3,3, 3]; % Quadratic B-spline in x-direction with 4 control points
knots_y = [-3,-3 -3, -1,1, 3,3, 3]; % Quadratic B-spline in y-direction with 4 control points

% Create the B-spline basis functions for x and y
n_x = length(knots_x) - px - 1;
n_y = length(knots_y) - py - 1;

x_int = [-1,1];
y_int = [-1,1];

t_x = linspace(-1, 1, 6000);  % Evaluation points in x
t_y = linspace(-1, 1, 6000);  % Evaluation points in y


B_x = zeros(length(t_x), n_x);
B_y = zeros(length(t_y), n_y);

for i = 1:n_x
    B_x(:, i) = utils().bspline_basis(i, px, knots_x, t_x);
end

for i = 1:n_y
    B_y(:, i) = utils().bspline_basis(i, py, knots_y, t_y);
end

% Create the design matrix for 2D B-spline
X = zeros(length(x), n_x * n_y);

for i = 1:length(x)
    for j = 1:n_x
        for k = 1:n_y
            X(i, (j - 1) * n_y + k) = utils().bspline_basis(j, px, knots_x, x(i)) * utils().bspline_basis(k, py, knots_y, y(i));
        end
    end
end

% Fit the 2D B-spline using least squares
C = X \ z';

% Evaluate the fitted 2D B-spline surface
fitted_surface = zeros(length(t_x), length(t_y));

for i = 1:length(t_x)
    for j = 1:length(t_y)
        for k = 1:n_x
            for l = 1:n_y
                fitted_surface(i, j) = fitted_surface(i, j) + C((k - 1) * n_y + l) * B_x(i, k) * B_y(j, l);
            end
        end
    end
end
[x_grid, y_grid] = meshgrid(t_x, t_y);
real=F(y_grid,x_grid);
Einf =norm(real(:)-fitted_surface(:),'Inf');
E2 =norm(real(:)-fitted_surface(:),2)/sqrt(numel(real));
disp(['Einf in the interval: ', num2str(Einf, '%.3e')]);
disp(['E2 in the interval: ', num2str(E2, '%.3e')]);