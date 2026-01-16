function u = utils()
% COMMON_PLOT_UTILS Provides common functions for visualizations
%
% This function returns a structure with common functions used
% in different visualization techniques to avoid code duplication.
%
% Returns:
%   utils: Structure with the following functions:
%     - apply_format: Applies a common format to figures.
%     - find_extraordinary_vertices: Finds vertices with valence different from 6.
%     - generate_figures: Generates figures for the initial and final states of a mesh.
%     - my_saveas: Saves figures in EPS format.
%     - make_draw_ring: Creates a TikZ figure of a ring with vertices inside and outside.
%     - make_figure_ring: Creates a MATLAB figure of a ring with vertices inside and outside.
%     - new_figure: Creates a new figure with a specified position.
%     - plotMesh: Plots a mesh given vertices and faces.
%     - make_draw: Creates a TikZ figure of a mesh with labeled vertices.
%     - basic_limit_function_figures_fast: Generates figures for the basic
%     limit function. (faster but slightly poorer visualiztion)
%     - basic_limit_function_figures: Generates figures for the basic limit function.
%     - minimum_bandwidth: Calculates the minimum bandwidth for a given mesh.
%     - measure_errors: Measures errors in a specified interval.
%     - locfit: Local regression and likelihood fitting function.
%     - predict: Predicts values using a fitted model.
%     - distance_matrix: Accumulate sum of squares of coordinate differences
%     - bspline_basis: Generates B-spline basis functions.
%     - get_noise_reduction_coefficient: Compute the noise reduction factor


u.apply_format = @apply_format;
u.find_extraordinary_vertices = @find_extraordinary_vertices;
u.generate_figures = @generate_figures;
u.my_saveas = @my_saveas;
u.make_draw_ring = @make_draw_ring;
u.make_figure_ring = @make_figure_ring;
u.new_figure = @new_figure;
u.plotMesh = @plotMesh;
u.make_draw = @make_draw;
u.basic_limit_function_figures_fast = @basic_limit_function_figures_fast;
u.basic_limit_function_figures = @basic_limit_function_figures;
u.minimum_bandwidth = @minimum_bandwidth;
u.measure_errors = @measure_errors;
u.locfit = @locfit;
u.predict = @predict;
u.distance_matrix = @distance_matrix;
u.bspline_basis = @bspline_basis;
u.get_noise_reduction_coefficient = @get_noise_reduction_coefficient
end

function apply_format(h,grid_on,num_grafica)
if nargin < 3
    num_grafica = 0;
end

figure(h);
colormap default
axis tight;
axis square;
axis off;
view(3);

if grid_on
    h.Children(1).Children(1).LineWidth = 2;
else
    h.Children(1).Children(1).LineStyle = 'none';
end
switch num_grafica
    case 1
        ll = light;
        ll.Position = [-1 -1 1];
    case 2
        ll = light;
        ll.Position = [-1 -1 1];
        axis([0.2 0.8 0.2 0.8])
    case 4
        axis([-2 2 -2 2])
        ll = light;
        ll.Position = [-1 -1 1];
    case 5
        axis([-1 1 -1 1])
        view(0,90)
    case 6
        ll = light;
        ll.Position = [-1 -1 1];
    case 7
        view(0,90)
        axis on;
        set(gca,'FontSize',20);
        % get axis limits
        xlim = get(gca,'XLim');
        ylim = get(gca,'YLim');
        % ceil them
        xlim = [floor(xlim(1))-1e-4, ceil(xlim(2))+1e-4];
        ylim = [floor(ylim(1))-1e-4, ceil(ylim(2))+1e-4];
        % set them back
        set(gca,'XLim',xlim);
        set(gca,'YLim',ylim);
        axis equal;
end
end

function extraordinary_vertices = find_extraordinary_vertices(faces)
% Find the list of vertices with valence different from 6
% Input:
%   faces - 3xN matrix with the faces of a triangulation
% Output:
%   extraordinary_vertices - list of vertices with valence different from 6

% Find the number of vertices
num_vertices = max(faces(:));

% Initialize valence array
valence = zeros(1, num_vertices);

% Calculate the valence of each vertex
for i = 1:size(faces, 2)
    for j = 1:3
        valence(faces(j, i)) = valence(faces(j, i)) + 1;
    end
end

% Find vertices with valence different from 6
extraordinary_vertices = find(valence ~= 6);
end

function generate_figures(vertices, faces,vertices_0,faces_0,gn,save_results,filename,N)

h = new_figure;
scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),50,'filled');
apply_format(h,true,gn);
if save_results
    my_saveas(h,[filename,'_',num2str(0),'iter_dots']); %#ok<*UNRCH>
end

h = new_figure;
plotMesh(vertices_0, faces_0);
apply_format(h,true,gn);
if save_results
    my_saveas(h,[filename,'_',num2str(0),'iter']); %#ok<*UNRCH>
end

h = new_figure;
plotMesh(vertices, faces);
apply_format(h,false, gn);
if save_results
    my_saveas(h,[filename,'_',num2str(N),'iter']);
end
hold on;
scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),50,'filled');
if save_results
    my_saveas(h,[filename,'_',num2str(N),'iter_dots']);
end
end

function my_saveas(h,nombre,pos)
if nargin<3
    pos = [1           1        1440         900];
end
h.Position=pos;
figure(h);
drawnow;
saveas(h,nombre,'epsc');
% exportgraphics(h,[nombre,'.eps']);
end

function make_draw_ring(filename,vertices,faces,center, radius, labels)
if nargin < 6
    labels = {};
end
mini = min(vertices(1:2,:),[],2);
width = max(vertices(2,:)) - min(vertices(2,:));
width = 0.5*max(1e-8,width);

% get the vertex indices that are inside the ring
distances = sqrt(sum((vertices(1:2,:) - center(1:2)).^2, 1));
inside = distances < max(radius)-1e-15;

fileID = fopen([filename,'.txt'],'w');
fprintf(fileID,'\\begin{tikzpicture}\n');

vertices = vertices(1:2,:);
vertices = (vertices-mini)/width;
vertex_in = vertices(:,inside);
vertex_out = vertices(:,~inside);
center = center(1:2);
center = (center-mini)/width;
radius = radius / width;

% draw a green circle.
fprintf(fileID,'\\draw[green] (%f,%f) circle (%f);\n',center(1),center(2),radius(1));
if length(radius) > 1
    % draw a blue circle.
    fprintf(fileID,'\\draw[blue] (%f,%f) circle (%f);\n',center(1),center(2),radius(2));
end

order = [1:size(faces,1),1];
edges = '\\draw[line width=0.05mm,   black] (%f,%f) -- (%f,%f);\n';
for i=1:size(faces,2)
    for j=1:size(faces,1)
        fprintf(fileID,edges,...
            vertices(1,faces(order(j),i)),...
            vertices(2,faces(order(j),i)),...
            vertices(1,faces(order(j+1),i)),...
            vertices(2,faces(order(j+1),i)));
    end
end

points_red = '\\filldraw[red] (%f,%f) circle (0.1pt);\n';
for i=1:size(vertex_in,2)
    fprintf(fileID,points_red,...
        vertex_in(1,i),vertex_in(2,i));
end

if ~isempty(labels)
    points_labels = '\\node[anchor=south west, font=\\tiny] at (%.2f,%.2f) {%s};\n';
    % sorting vertices: higher Y first, higher X first.
    [~, idx] = sortrows(vertex_in', [-2, -1]); % Sort by Y descending, then X descending
    vertex_in_order = vertex_in(:, idx); % Reorder vertices based on sorted indices
    for i=1:size(vertex_in_order,2)
        fprintf(fileID,points_labels,...
            vertex_in_order(1,i)-0.1,vertex_in_order(2,i)-0.1,labels{i});
    end
end

points_black = '\\filldraw[black] (%f,%f) circle (0.1pt);\n';
for i=1:size(vertex_out,2)
    fprintf(fileID,points_black,vertex_out(1,i),vertex_out(2,i));
end

% draw a red cross at the center
fprintf(fileID,'\\draw[red,line width=0.1mm] (%f,%f) -- (%f,%f);\n',center(1)-0.02,center(2)-0.02,center(1)+0.02,center(2)+0.02);
fprintf(fileID,'\\draw[red,line width=0.1mm] (%f,%f) -- (%f,%f);\n',center(1)-0.02,center(2)+0.02,center(1)+0.02,center(2)-0.02);

fprintf(fileID,'\\end{tikzpicture}');
fclose(fileID);
end

function make_figure_ring(vertices, faces, center, radius,labels)
figure;
hold on;
plotMesh(vertices,faces,3);
view([0,90]);
axis equal;
axis off;

% get the vertex indices that are inside the ring
distances = sqrt(sum((vertices(1:2,:) - center(1:2)).^2, 1));
inside = distances < max(radius)-1e-15;
vertex_in = vertices(:,inside);
vertex_out = vertices(:,~inside);


% show dots at the vertices
plot(vertex_in(1,:),vertex_in(2,:),'r.','MarkerSize',20);
plot(vertex_out(1,:),vertex_out(2,:),'k.','MarkerSize',20);

% sorting vertices: higher Y first, higher X first.
[~, idx] = sortrows(vertex_in', [-2, -1]); % Sort by Y descending, then X descending
vertex_in_order = vertex_in(:, idx); % Reorder vertices based on sorted indices

% labels for the inner vertices, shown in the top right of the vertices
for i = 1:length(labels)
    text(vertex_in_order(1,i), vertex_in_order(2,i), labels{i}, 'Color', 'black', 'FontSize', 50, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

% mark the center with a red cross
plot(center(1), center(2), 'rx', 'MarkerSize', 15, 'LineWidth', 5);

% show green circle
theta = linspace(0, 2*pi, 100);
x = center(1) + radius(1) * cos(theta);
y = center(2) + radius(1) * sin(theta);
plot(x, y, 'g-', 'LineWidth', 2);

if length(radius) > 1
    % show blue circle
    theta = linspace(0, 2*pi, 100);
    x = center(1) + radius(2) * cos(theta);
    y = center(2) + radius(2) * sin(theta);
    plot(x, y, 'b-', 'LineWidth', 2);
end
end

function h = new_figure
h = figure;
% asd=get(groot,'MonitorPositions');
% asd = asd(end,:);
h.Position=[1           1        1440         900];
end

function plotMesh(vertices, faces,op)
vertices=double(vertices);
hold on;
if nargin < 3
    op = 1; % default operation
elseif numel(op) > 1
    X = vertices;
    Y = faces;
    Z = op;
    op = 4;
end
switch op
    case 1
        s=trisurf(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        %             light('Position',[1 1 1])
        %             s.EdgeColor = 'none';
        s.EdgeColor = [0,0,0];
        s.FaceColor = 'interp';
        s.LineWidth = 2;
        %             colormap gray(1);
        colormap default
    case 2
        trimesh(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        colormap default
        %     colormap gray(1);
    case 3 % transparent faces
        s = trisurf(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        s.EdgeColor = [0,0,0];
        s.LineWidth = 2;
        s.FaceColor = 'none';
    case 4
        s=surf(X, Y, Z);
        s.EdgeColor = 'none';
        s.FaceColor = 'interp';
        s.LineWidth = 2;
        colormap default
end
axis tight;
axis square;
axis off;
view(3);
end

function make_draw(filename,vertices,faces,labels)
% vertices contains the 2D coordinates
% faces contains the indices of the vertices that are connected
% labels is a cell containing string labels for each vertex
mini = min(vertices,[],"all");
maxi = max(vertices,[],"all");
maxi = max(mini+1e-8,maxi);
vertices = 10*(vertices-mini)/(maxi-mini);
edges = '\\draw[line width=0.3mm,   black] (%.2f,%.2f) -- (%.2f,%.2f);\n';
points = '\\filldraw[black] (%.2f,%.2f) circle (1.5pt);\n';
% point labels where no point or edge is drawn
points_labels = '\\node[anchor=south west] at (%.2f,%.2f) {%s};\n';

fileID = fopen([filename,'.txt'],'w');
fprintf(fileID,'\\begin{tikzpicture}\n');
order = [1:size(faces,1),1];
for i=1:size(faces,2)
    for j=1:size(faces,1)
        fprintf(fileID,edges,...
            vertices(1,faces(order(j),i)),...
            vertices(2,faces(order(j),i)),...
            vertices(1,faces(order(j+1),i)),...
            vertices(2,faces(order(j+1),i)));
    end
end
for i=1:size(vertices,2)
    fprintf(fileID,points,...
        vertices(1,i),vertices(2,i));
    fprintf(fileID,points_labels,...
        vertices(1,i),vertices(2,i),labels{i});
end
fprintf(fileID,'\\end{tikzpicture}');
fclose(fileID);
end

function basic_limit_function_figures_fast(W, L, basis, filename, save_results)
A = double(compute_mask(W,L,L,basis,false));

N = ceil(7*L);
f = zeros(2*N+1,2*N+1);
f(N+1,N+1)=1;
[X,Y] = meshgrid(-N:N,-N:N);
% apply basis
points = transpose([X(:),Y(:)]);
basis = transpose(basis);
points = basis*points;
X = reshape(points(1,:),size(X));
Y = reshape(points(2,:),size(Y));

it = 5;
f = apply_mask(f, A, it);
X = apply_mask(X, A, it);
Y = apply_mask(Y, A, it);

bool = f(:) > 0;
f(~bool) = NaN;

h = utils().new_figure();
utils().plotMesh(X,Y,f);
utils().apply_format(h,true,6);
if save_results
    utils().my_saveas(h,filename);
end

h = utils().new_figure();
utils().plotMesh(X,Y,f);
utils().apply_format(h,false,7);
if save_results
    utils().my_saveas(h,[filename,'_support']);
end
end

function basic_limit_function_figures(vertices, faces, L, weight, filename, save_results)
% find the [0;0;0] vertex
[~, idx] = find(vertices(1,:)==0 & vertices(2,:)==0 & vertices(3,:)==0);
vertices(3,idx) = 1;

vertices_0 = vertices;
faces_0 = faces;

N=5;
for i=1:N
    [vertices, faces] = RegressionSubdivision(vertices, faces, L/2^(i-1), weight,'functional');
end


h = utils().new_figure();
utils().plotMesh(vertices_0, faces_0);
utils().apply_format(h,true,6);
if save_results
    utils().my_saveas(h,[filename,'_',num2str(0),'iter']); %#ok<*UNRCH>
end

bool = vertices(3,:) > 0;
vertices(3,~bool) = NaN;

h = utils().new_figure();
utils().plotMesh(vertices, faces);
utils().apply_format(h,false,6);
if save_results
    utils().my_saveas(h,[filename,'_',num2str(N),'iter']);
end
save([filename,'_',num2str(N),'iter'],'vertices','faces')


h = utils().new_figure();
utils().plotMesh(vertices, faces);
utils().apply_format(h,false,7);
if save_results
    utils().my_saveas(h,[filename,'_support']);
end
end

function L = minimum_bandwidth(vertices, faces)
L = 0;
for i=1:size(faces,2)
    v1 = vertices(:,faces(1,i));
    v2 = vertices(:,faces(2,i));
    v3 = vertices(:,faces(3,i));
    dist12 = norm(v1-v2);
    dist13 = norm(v1-v3);
    dist23 = norm(v2-v3);
    L = max(L,max([dist12,dist13,dist23]));
end
end

function [Einf,E2] = measure_errors(vertices, F)
x_int = [-1,1];
y_int = [-1,1];
bool_int = (vertices(1,:)>=x_int(1)) & (vertices(1,:)<=x_int(2)) & (vertices(2,:)>=y_int(1)) & (vertices(2,:)<=y_int(2));
Einf = norm(vertices(3,bool_int)-F(vertices(1,bool_int),vertices(2,bool_int)),'Inf');
E2 = norm(vertices(3,bool_int)-F(vertices(1,bool_int),vertices(2,bool_int)),2)/sqrt(length(vertices(3,bool_int)));
disp(['Einf in the interval: ', num2str(Einf, '%.3e')]);
disp(['E2 in the interval: ', num2str(E2, '%.3e')]);
end

function fit=locfit(varargin)

% Smoothing noisy data using Local Regression and Likelihood.
%
% arguments still to add: dc maxit
%
%  Usage: fit = locfit(x,y)   % local regression fit of x and y.
%         fit = locfit(x)     % density estimation of x.
%
%  Smoothing with locfit is a two-step procedure. The locfit()
%  function evaluates the local regression smooth at a set of points
%  (can be specified through an evaluation structure). Then, use
%  the predict() function to interpolate this fit to other points.
%
%  Additional arguments to locfit() are specified as 'name',value pairs, e.g.:
%  locfit( x, 'alpha',[0.7,1.5] , 'family','rate' , 'ev','grid' , 'mg',100 );
%
%
%  Data-related inputs:
%
%    x is a vector or matrix of the independent (or predictor) variables.
%      Rows of x represent subjects, columns represent variables.
%      Generally, local regression would be used with 1-4 independent
%      variables. In higher dimensions, the curse-of-dimensionality,
%      as well as the difficulty of visualizing higher dimensional
%      surfaces, may limit usefulness.
%
%    y is the column vector of the dependent (or response) variable.
%      For density families, 'y' is omitted.
% NOTE: x and y are the first two arguments. All other arguments require
%        the 'name',value notation.
%
%    'weights' Prior weights for observations (reciprocal of variance, or
%           sample size).
%    'cens' Censoring indicators for hazard rate or censored regression.
%           The coding is '1' (or 'TRUE') for a censored observation, and
%           '0' (or 'FALSE') for uncensored observations.
%    'base' Baseline parameter estimate. If a baseline is provided,
%           the local regression model is fitted as
%                        Y_i = b_i + m(x_i) + epsilon_i,
%           with Locfit estimating the m(x) term. For regression models,
%           this effectively subtracts b_i from Y_i. The advantage of the
%           'base' formulation is that it extends to likelihood
%           regression models.
%    'scale' A scale to apply to each variable. This is especially
%           important for multivariate fitting, where variables may be
%           measured in non-comparable units. It is also used to specify
%           the frequency for variables with the 'a' (angular) style.
%     'sty' Character string (length d) of styles for each predictor variable.
%           n denotes `normal'; a denotes angular (or periodic); l and r
%           denotes one-sided left and right; c is conditionally parametric.
%
%
%  Smoothing Parameters and Bandwidths:
%  The bandwidth (or more accurately, half-width) of the smoothing window
%  controls the amount of smoothing. Locfit allows specification of constant
%  (fixed), nearest neighbor, certain locally adaptive variable bandwidths,
%  and combinations of these. Also related to the smoothing parameter
%  are the local polynmial degree and weight function.
%
%    'nn' 'Nearest neighbor' smoothing parameter. Specifying 'nn',0.5
%         means that the width of each smoothing neighborhood is chosen
%         to cover 50% of the data.
%
%     'h' A constant (or fixed) bandwidth parameter. For example, 'h',2
%         means that the smoothing windows have constant half-width
%         (or radius) 2. Note that h is applied after scaling.
%
%   'pen' penalty parameter for adaptive smoothing. Needs to be used
%         with care.
%
%  'alpha' The old way of specifying smoothing parameters, as used in
%         my book. alpha is equivalent to the vector [nn,h,pen].
%         If multiple componenents are non-zero, the largest corresponding
%         bandwidth is used. The default (if none of alpha,nn,h,pen
%         are provided) is [0.7 0 0].
%
%   'deg' Degree of local polynomial. Default: 2 (local quadratic).
%         Degrees 0 to 3 are supported by almost all parts of the
%         Locfit code. Higher degrees may work in some cases.
%
%  'kern' Weight function, default = 'tcub'. Other choices are
%         'rect', 'trwt', 'tria', 'epan', 'bisq' and 'gauss'.
%         Choices may be restricted when derivatives are
%         required; e.g. for confidence bands and some bandwidth
%         selectors.
%
%    'kt' Kernel type, 'sph' (default); 'prod'. In multivariate
%         problems, 'prod' uses a simplified product model which
%         speeds up computations.
%
%  'acri' Criterion for adaptive bandwidth selection.
%
%
%  Derivative Estimation.
%  Generally I recommend caution when using derivative estimation
%  (and especially higher order derivative estimation) -- can you
%  really estimate derivatives from noisy data? Any derivative
%  estimate is inherently more dependent on an assumed smoothness
%  (expressed through the bandwidth) than the data. Warnings aside...
%
%  'deriv' Derivative estimation. 'deriv',1 specifies the first derivative
%         (or more correctly, an estimate of the local slope is returned.
%         'deriv',[1 1] specifies the second derivative. For bivariate fits
%         'deriv',2 specifies the first partial derivative wrt x2.
%         'deriv',[1 2] is mixed second-order derivative.
%
%  Fitting family.
%  'family' is used to specify the local likelihood family.
%         Regression-type families are 'gaussian', 'binomial',
%           'poisson', 'gamma' and 'geom'. If the family is preceded
%           by a q (e.g. 'qgauss', or 'qpois') then quasi-likelihood is
%           used; in particular, a dispersion estimate is computed.
%           Preceding by an 'r' makes an attempt at robust (outlier-resistant)
%           estimation. Combining q and r (e.g. 'family','qrpois') may
%           work, if you're lucky.
%         Density estimation-type families are 'dens', 'rate' and 'hazard'
%           (hazard or failure rate). Note that `dens' scales the output
%           to be a statistical density estimate (i.e. scaled to integrate
%           to 1). 'rate' estimates the rate or intensity function (events
%           per unit time, or events per unit area), which may be called
%           density in some fields.
%         The default family is 'qgauss' if a response (y argument) has been
%         provided, and 'dens' if no response is given.
%    'link' Link function for local likelihood fitting. Depending on the
%           family, choices may be 'ident', 'log', 'logit',
%           'inverse', 'sqrt' and 'arcsin'.
%
%  Evaluation structures.
%    By default, locfit chooses a set of points, depending on the data
%    and smoothing parameters, to evaluate at. This is controlled by
%    the evaluation structure.
%      'ev' Specify the evaluation structure. Default is 'tree'.
%           Other choices include 'phull' (triangulation), 'grid' (a grid
%           of points), 'data' (each data point), 'crossval' (data,
%           but use leave-one-out cross validation), 'none' (no evaluation
%           points, effectively producing the global parametric fit).
%           Alternatively, a vector/matrix of evaluation points may be
%           provided.
%           (kd trees not currently supported in mlocfit)
%     'll' and 'ur' -- row vectors specifying the upper and lower limits
%           for the bounding box used by the evaluation structure.
%           They default to the data range.
%     'mg' For the 'grid' evaluation structure, 'mg' specifies the
%           number of points on each margin. Default 10. Can be either a
%           single number or vector.
%    'cut' Refinement parameter for adaptive partitions. Default 0.8;
%           smaller values result in more refined partitions.
%    'maxk' Controls space assignment for evaluation structures. For the
%           adaptive evaluation structures, it is impossible to be sure
%           in advance how many vertices will be generated. If you get
%           warnings about `Insufficient vertex space', Locfit's default
%           assigment can be increased by increasing 'maxk'. The default
%           is 'maxk','100'.
%
%    'xlim' For density estimation, Locfit allows the density to be
%           supported on a bounded interval (or rectangle, in more than
%           one dimension). The format should be [ll;ul] (ie, matrix with
%           two rows, d columns) where ll is the lower left corner of
%           the rectangle, and ur is the upper right corner.
%           One-sided bounds, such as [0,infty), are not supported, but can be
%           effectively specified by specifying a very large upper
%           bound.
%
%      'module' either 'name' or {'name','/path/to/module',parameters}.
%
%  Density Estimation
%      'renorm',1  will attempt to renormalize the local likelihood
%           density estimate so that it integrates to 1. The llde
%           (specified by 'family','dens') is scaled to estimate the
%           density, but since the estimation is pointwise, there is
%           no guarantee that the resulting density integrates exactly
%           to 1. Renormalization attempts to achieve this.
%
%  The output of locfit() is a Matlab structure:
%
% fit.data.x (n*d)
% fit.data.y (n*1)
% fit.data.weights (n*1 or 1*1)
% fit.data.censor (n*1 or 1*1)
% fit.data.baseline (n*1 or 1*1)
% fit.data.style (string length d)
% fit.data.scales (1*d)
% fit.data.xlim (2*d)
%
% fit.evaluation_structure.type (string)
% fit.evaluation_structure.module.name (string)
% fit.evaluation_structure.module.directory (string)
% fit.evaluation_structure.module.parameters (string)
% fit.evaluation_structure.lower_left (numeric 1*d)
% fit.evaluation_structure.upper_right (numeric 1*d)
% fit.evaluation_structure.grid (numeric 1*d)
% fit.evaluation_structure.cut (numeric 1*d)
% fit.evaluation_structure.maxk
% fit.evaluation_structure.derivative
%
% fit.smoothing_parameters.alpha = (nn h pen) vector
% fit.smoothing_parameters.adaptive_criterion (string)
% fit.smoothing_parameters.degree (numeric)
% fit.smoothing_parameters.family (string)
% fit.smoothing_parameters.link (string)
% fit.smoothing_parameters.kernel (string)
% fit.smoothing_parameters.kernel_type (string)
% fit.smoothing_parameters.deren
% fit.smoothing_parameters.deit
% fit.smoothing_parameters.demint
% fit.smoothing_parameters.debug
%
% fit.fit_points.evaluation_points (d*nv matrix)
% fit.fit_points.fitted_values (matrix, nv rows, many columns)
% fit.fit_points.evaluation_vectors.cell
% fit.fit_points.evaluation_vectors.splitvar
% fit.fit_points.evaluation_vectors.lo
% fit.fit_points.evaluation_vectors.hi
% fit.fit_points.fit_limits (d*2 matrix)
% fit.fit_points.family_link (numeric values)
% fit.fit_points.kappa (likelihood, degrees of freedom, etc)
%
% fit.parametric_component
%
%
%  The OLD format:
%
%    fit{1} = data.
%    fit{2} = evaluation structure.
%    fit{3} = smoothing parameter structure.
%    fit{4}{1} = fit points matrix.
%    fit{4}{2} = matrix of fitted values etc.
%           Note that these are not back-transformed, and may have the
%           parametric component removed.
%           (exact content varies according to module).
%    fit{4}{3} = various details of the evaluation points.
%    fit{4}{4} = fit limits.
%    fit{4}{5} = family,link.
%    fit{5} = parametric component values.
%



% Minimal input validation
if nargin < 1
    error( 'At least one input argument required' );
end

xdata = double(varargin{1});
d = size(xdata,2);
n = size(xdata,1);
if ((nargin>1) && (~ischar(varargin{2})))
    ydata = double(varargin{2});
    if (any(size(ydata) ~= [n 1]));
        error('y must be n*1 column vector');
    end
    family = 'qgauss';
    na = 3;
else
    ydata = 0;
    family = 'density';
    na = 2;
end
if mod(nargin-na,2)==0
    error( 'All arguments other than x, y must be name,value pairs' );
end


wdata = ones(n,1);
cdata = 0;
base  = 0;
style = 'n';
scale = 1;
xl = zeros(2,d);

alpha = [0 0 0];
deg = 2;
link = 'default';
acri = 'none';
kern = 'tcub';
kt = 'sph';
deren = 0;
deit  = 'default';
demint= 20;
debug = 0;

ev = 'tree';
ll = zeros(1,d);
ur = zeros(1,d);
mg = 10;
maxk = 100;
deriv=0;
cut = 0.8;
mdl = struct('name','std', 'directory','', 'parameters',0 );

while na < length(varargin)
    inc = 0;
    if (varargin{na}=='y')
        ydata = double(varargin{na+1});
        family = 'qgauss';
        inc = 2;
        if (any(size(ydata) ~= [n 1]))
            error('y must be n*1 column vector');
        end
    end
    if (strcmp(varargin{na},'weights'))
        wdata = double(varargin{na+1});
        inc = 2;
        if (any(size(wdata) ~= [n 1]))
            error('weights must be n*1 column vector');
        end
    end
    if (strcmp(varargin{na},'cens'))
        cdata = double(varargin{na+1});
        inc = 2;
        if (any(size(cdata) ~= [n 1]))
            error('cens must be n*1 column vector');
        end
    end
    if (strcmp(varargin{na},'base')) % numeric vector, n*1 or 1*1.
        base = double(varargin{na+1});
        if (length(base)==1)
            base = base*ones(n,1);
        end
        inc = 2;
    end
    if (strcmp(varargin{na},'style')) % character string of length d.
        style = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'scale')) % row vector, length 1 or d.
        scale = varargin{na+1};
        if (scale==0)
            scale = zeros(1,d);
            for i=1:d
                scale(i) = sqrt(var(xdata(:,i)));
            end
        end
        inc = 2;
    end
    if (strcmp(varargin{na},'xlim')) % 2*d numeric matrix.
        xl = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'alpha')) % row vector of length 1, 2 or 3.
        alpha = [varargin{na+1} 0 0 0];
        alpha = alpha(1:3);
        inc = 2;
    end
    if (strcmp(varargin{na},'nn')) % scalar
        alpha(1) = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'h')) % scalar
        alpha(2) = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'pen')) % scalar
        alpha(3) = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'acri')) % string
        acri = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'deg')) % positive integer.
        deg = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'family')) % character string.
        family = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'link')) % character string.
        link = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'kern')) % character string.
        kern = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'kt')) % character string.
        kt = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'ev')) % char. string, or matrix with d columns.
        ev = varargin{na+1};
        if (isnumeric(ev))
            ev = ev';
        end
        inc = 2;
    end
    if (strcmp(varargin{na},'ll')) % row vector of length d.
        ll = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'ur')) % row vector of length d.
        ur = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'mg')) % row vector of length d.
        mg = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'cut')) % positive scalar.
        cut = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'module')) % string.
        mdl = struct('name',varargin{na+1}, 'directory','', 'parameters',0 );
        inc = 2;
    end
    if (strcmp(varargin{na},'maxk')) % positive integer.
        maxk = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'deriv')) % numeric row vector, up to deg elements.
        deriv = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'renorm')) % density renormalization.
        deren = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'itype')) % density - integration type.
        deit = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'mint')) % density - # of integration points.
        demint = varargin{na+1};
        inc = 2;
    end
    if (strcmp(varargin{na},'debug')) % debug level.
        debug = varargin{na+1};
        inc = 2;
    end
    if (inc==0)
        disp(varargin{na});
        error('Unknown Input Argument.');
    end
    na=na+inc;
end


fit.data.x = xdata;
fit.data.y = ydata;
fit.data.weights = wdata;
fit.data.censor = cdata;
fit.data.baseline = base;
fit.data.style = style;
fit.data.scales = scale;
fit.data.xlim = xl;

fit.evaluation_structure.type = ev;
fit.evaluation_structure.module = mdl;
fit.evaluation_structure.lower_left = ll;
fit.evaluation_structure.upper_right = ur;
fit.evaluation_structure.grid = mg;
fit.evaluation_structure.cut = cut;
fit.evaluation_structure.maxk = maxk;
fit.evaluation_structure.derivative = deriv;

if (alpha==0)
    alpha = [0.7 0 0];
end

fit.smoothing_parameters.alpha = alpha;
fit.smoothing_parameters.adaptive_criterion = acri;
fit.smoothing_parameters.degree = deg;
fit.smoothing_parameters.family = family;
fit.smoothing_parameters.link = link;
fit.smoothing_parameters.kernel = kern;
fit.smoothing_parameters.kernel_type = kt;
fit.smoothing_parameters.deren = deren;
fit.smoothing_parameters.deit = deit;
fit.smoothing_parameters.demint = demint;
fit.smoothing_parameters.debug = debug;

[fpc pcomp] = mexlf(fit.data,fit.evaluation_structure,fit.smoothing_parameters);
fit.fit_points = fpc;
fit.parametric_component = pcomp;

return
end

function [y, se] = predict(varargin)

% Interpolate a fit produced by locfit().
%
% predict(fit)    produces the fitted values at locfit's selected points.
% predict(fit,x)  interpolates the fits to points specified by x.
%
% Input arguments:
%   fit   The locfit() fit.
%   x     Points to interpolate at. May be a matrix with d columns,
%         or cell with d components (each a vector). In the former
%         case, a fitted value is computed for each row of x.
%         In the latter, the components of x are interpreted as
%         grid margins.
%         Can also specify 'data' (evaluate at data points);
%         or 'fitp' (extract the fitted points).
%  'band',value
%         Type of standard errors to compute. Default is 'band','n', for none.
%         Other choices are 'band','g' (use a global s to estimate the resiudal
%         standard deviation, so standard errors are s*||l(x)||);
%         'band','l' (use a local s(x), so std. errors are s(x)*||l(x)||);
%         'band','p' (prediction errors, so s*sqrt(1+||l(x)||^2).
%  'direct'
%         Compute the local fit directly (rather than using local
%         regression, at each point specified by the x argument.
%  'kappa',vector
%         Vector of constants for simultaneous confidence bands,
%         computed by the kappa0() function.
%  'level',value
%         Coverage probability for confidence intervals and bands.
%         Default is 0.95.
%
%  Output is a vector of fitted values (if 'band','n'), or a cell
%  with fitted value, standard error vectors, and matrix of lower
%  and upper confidence limits.
%
%  Note that for local likelihood fits, back-transformation is
%  not performed, so that (e.g.) for Poisson regression with the
%  log-link, the output estimates the log-mean, and its standard errors.
%  Likewise, for density estimation, the output is log(density).
%
%  Author: Catherine Loader.

if (nargin<1)
    error('predict requires fit argument');
end

fit = varargin{1};

if (nargin==1)
    x = 'fitp';
else
    x = varargin{2};
end

band = 'n';
what = 'coef';
rest = 'none';
dir  = 0;
level = 0.95;
d = size(fit.data.x,2);
kap = [zeros(1,d) 1];

na = 3;
while na <= nargin
    inc = 0;
    if strcmp(varargin{na},'band')
        band = varargin{na+1};
        inc = 2;
    end
    if strcmp(varargin{na},'what')
        what = varargin{na+1};
        inc = 2;
    end
    if strcmp(varargin{na},'restyp')
        rest = varargin{na+1};
        inc = 2;
    end
    if strcmp(varargin{na},'direct')
        dir = 1;
        inc = 1;
    end
    if strcmp(varargin{na},'kappa')
        kap = varargin{na+1};
        inc = 2;
    end
    if strcmp(varargin{na},'level')
        level = varargin{na+1};
        inc = 2;
    end
    if (inc == 0)
        disp(varargin{na});
        error('Unknown argument');
    end
    na = na+inc;
end

[y se cb] = mexpp(x,fit,band,what,rest,dir,kap,level);
if (band=='n')
    y = y;
else
    y = {y se cb};
end

return;
end

function DM = distance_matrix(dsites,ctrs)
[M,s] = size(dsites); [N,s] = size(ctrs);
DM = zeros(M,N);
% Accumulate sum of squares of coordinate differences
% The ndgrid command produces two MxN matrices:
%'/„ . dr, consisting of N identical columns (each containing
% the d-th coordinate of the M data sites)
% cc, consisting of M identical rows (each containing
% the d-th coordinate of the N centers)
for d=1:s
    [dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);
end

function result = bspline_basis(i, p, knots, x)
if p == 0
    result = (knots(i) <= x).* (x < knots(i + 1));
else
    result = 0;
    if knots(i + p) ~= knots(i)
        result = result + (x - knots(i))./ (knots(i + p) - knots(i)).* bspline_basis(i, p - 1, knots, x);
    end
    if knots(i + p + 1) ~= knots(i + 1)
        result = result + (knots(i + p + 1) - x) ./ (knots(i + p + 1) - knots(i + 1)).* bspline_basis(i + 1, p - 1, knots, x);
    end
end
end

function theta = get_noise_reduction_coefficient(w,L,L_stencil,basis)

    A = compute_mask(w,L,L_stencil,basis,false);
    theta = 0;

    for i=1:2
        for j=1:2
            a = A(i:2:end,j:2:end);
            theta = max(theta,sum(a(:).^2));
        end
    end
end