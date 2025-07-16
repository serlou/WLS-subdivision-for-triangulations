% generates matlab figure

clear; close all;
load('datasets/dataset2.mat')
save_results = false;
wei = 'gaussian';
switch wei
    case 'hat'
        weight = @(x) 1-x;
    case 'gaussian'
        weight = @(x) exp(-(2.5*x).^2/2);
end
L = 2;

F=@(xx,yy) sin(xx).*cos(yy);
vertices = vertices_0;
vertices_0 = vertices;
faces_0 = faces;

N=5;
for i=1:N
    [vertices, faces] = RegressionSubdivision(vertices, faces, L/2^(i-1),weight,'functional');
end

h = utils().new_figure();
utils().plotMesh(vertices_0, faces_0);
utils().apply_format(h,true,4);
if save_results
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',wei,'_',num2str(0),'iter']); %#ok<*UNRCH> 
end
h = utils().new_figure();
utils().plotMesh(vertices, faces);
utils().apply_format(h,false,4);
if save_results
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',wei,'_',num2str(N),'iter']);
end
hold on;
scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),25,'filled');
if save_results
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',wei,'_',num2str(N),'iter_dots']);
    save(['graphics/grafica4_sin_cos_',wei,'_',num2str(N),'iter'],'vertices','faces')
end

%%% To calculate the errors
utils().measure_errors(vertices, F);