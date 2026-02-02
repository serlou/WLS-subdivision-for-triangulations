clear;
close all;
load('datasets/dataset2.mat');
save_results = false;

F = @(x,y) sin(x).*cos(y);

%Table 1. Line 4
% W = @(e,r) exp(-(e*r).^2/2); ep = 2.5;

%Table 1. Line 10
W = @(e,r) exp(-(e*r).^2/2); ep = 2.5/2;

neval = 500;
[X,Y] = meshgrid([-3:3],[-3:3]);
x(:,1)=reshape(X,49,1);
x(:,2)=reshape(Y,49,1);
ctrs = x;
N=2;
M=49;
xx=vertices_0(1:2,:)';
dsites=xx;
z=vertices_0(3,:)';
DM_data = utils().distance_matrix(dsites,ctrs);
CM = W(ep,DM_data);   % Collocation matrix
rhs=z;
grid = linspace(-2,2,neval); [xe,ye] = meshgrid(grid);
epoints = [xe(:) ye(:)];
% Compute distance matrix between evaluation points and centers
DM_eval = utils().distance_matrix(epoints,ctrs);
EM = W(ep,DM_eval);   % Evaluation matrix
% Compute W least squares approximation
Pf = EM * (CM\rhs);
h1=utils().new_figure();
scatter3(xx(:,1),xx(:,2),rhs,'filled');
hold on
surf(xe,ye,reshape(Pf,neval,neval))
utils().apply_format(h1,false,4);
if save_results
    utils().my_saveas(h1,'graphics/rbfgaussep2.eps');
end

%%% To calculate the errors
load('datasets/dataset2_after_scheme.mat');
epoints = [vertices(1,:)' vertices(2,:)'];
% Compute distance matrix between evaluation points and centers
DM_eval = utils().distance_matrix(epoints,ctrs);
EM = W(ep,DM_eval);   % Evaluation matrix
% Add columns for linear reproduction
%PM = [ones(neval^2,1) epoints]; EM = [EM PM];
% Compute W least squares approximation
Pf = EM * (CM\rhs);
utils().measure_errors([vertices(1:2,:); Pf'], F);