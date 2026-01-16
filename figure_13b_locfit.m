clear;
close all;
load('datasets/dataset2.mat')
save_results = false;
xx=vertices_0(1:2,:)';
z=vertices_0(3,:)';
n=length(z);
%% Table 1. line 1
fit = utils().locfit(xx,z,'alpha',[0,1],'deg',1,'kern','tria');

%%Table 1. line 2
% fit = utils().locfit(xx,z,'alpha',[0,1],'deg',1,'kern','gauss');

%%Table 1. line 7
% fit = utils().locfit(xx,z,'alpha',[0,2],'deg',1,'kern','tria');

%%Table 1. line 8
% fit = utils().locfit(xx,z,'alpha',[0,2],'deg',1,'kern','gauss');

neval=100;
grid = linspace(-2,2,neval); [xe,ye] = meshgrid(grid);
Z1=utils().predict(fit,[xe(:) ye(:)]);
h1=utils().new_figure();
scatter3(xx(:,1),xx(:,2),z,'filled')
F=@(x,y)sin(x).*cos(y);
hold on
s1 = surf(xe,ye,reshape(Z1,neval,neval));
utils().apply_format(h1,false,4);
if save_results
    utils().my_saveas(h1,'graphics/llr1gaussr0alambda1.eps');
end

%%% To calculate the errors
load('datasets/dataset2_after_scheme.mat');
Pf=utils().predict(fit,[vertices(1,:)' vertices(2,:)']);
utils().measure_errors([vertices(1:2,:); Pf'], F);