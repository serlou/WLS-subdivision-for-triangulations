clear;
close all;
load('datasets/dataset2.mat');
save_results = false;
xx=vertices_0(1:2,:)';
z=vertices_0(3,:)';
F=@(x,y)sin(x).*cos(y);

%% Table 1. line 3
W = @(e,r) exp(-(e*r).^2/2); ep = 2.5;

%% Table 1. line 9
% W = @(e,r) exp(-(e*r).^2/2); ep = 2.5/2;

n=length(z);
neval=500;
grid = linspace(-2,2,neval); [xe,ye] = meshgrid(grid);
neval=length(xe);
epoints = [xe(:) ye(:)];
ctrs = xx;

DM_eval = utils().distance_matrix(epoints,ctrs);

EM = W(ep,DM_eval);
EM = EM./repmat(EM*ones(n,1),1,n); 
Pf = EM*z;

h1=utils().new_figure();
scatter3(xx(:,1),xx(:,2),z,'filled')
hold on
surf(xe,ye,reshape(Pf,neval,neval))
utils().apply_format(h1,false,4);
if save_results
    utils().my_saveas(h1,'graphics/shepard.eps');
end

%%% To calculate the errors
load('datasets/dataset2_after_scheme.mat');
epoints = [vertices(1,:)' vertices(2,:)'];
DM_eval = utils().distance_matrix(epoints,ctrs);
EM = W(ep,DM_eval);
EM = EM./repmat(EM*ones(n,1),1,n); % Shepard normalization °/0 Compute quasi-interpolant
Pf = EM*z;
utils().measure_errors([vertices(1:2,:); Pf'], F);