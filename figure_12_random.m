% generates matlab figure

clear; close all;
load('datasets/dataset2.mat')
save_results = true;
wei = 'hat';
switch wei
    case 'hat'
        weight = @(x) 1-x;
    case 'gaussian'
        weight = @(x) exp(-(2.5*x).^2/2);
end
L = 2;

rng(1);
vertices_0(3,:) = rand(size(vertices_0(3,:)));

vertices = vertices_0;
vertices_0 = vertices;
faces_0 = faces;

N=5;
for i=1:N
    [vertices, faces] = RegressionSubdivision(vertices, faces, L/2^(i-1),weight,'functional');
end

h1 = utils().new_figure();
utils().plotMesh(vertices_0, faces_0);
utils().apply_format(h1,true,0);
view(0,90)

extraordinary = utils().find_extraordinary_vertices(faces_0);

hold on;
scatter3(vertices_0(1,extraordinary), vertices_0(2,extraordinary), vertices_0(3,extraordinary)+1, 50, 'r', 'filled');
hold off;
axis on;
set(gca, 'FontSize', 30);
if save_results
    utils().my_saveas(h1,'graphics/grafica7_random_0iter');
end

h = utils().new_figure();
utils().plotMesh(vertices, faces);
utils().apply_format(h,false,0);
ll = light;
ll.Position = [0.2438 -0.1395 0.9476];
colormap(h, colormap(h1)); 
if save_results
    utils().my_saveas(h,['graphics/grafica7_random_',wei,'_',num2str(N),'iter']);
end
% hold on;
% scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),25,'filled');
% if save_results
%     utils().my_saveas(h,['graphics/grafica7_random_',wei,'_',num2str(N),'iter_dots']);
%     save(['graphics/grafica7_random_',wei,'_',num2str(N),'iter'],'vertices','faces')
% end