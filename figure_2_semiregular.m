% generates matlab figure

clear; close all;
load('datasets/dataset2.mat')
save_results = true;

vertices = vertices_0;
faces_0 = faces;

weight = @(x) 1;
L=2;

N=1;
for i=1:N
    [vertices, faces] = RegressionSubdivision(vertices, faces,L,weight,'functional');
end

h = utils().new_figure();
hold on;

extra = utils().find_extraordinary_vertices(faces_0);
scatter3(vertices_0(1,extra),vertices_0(2,extra),vertices_0(3,extra)+2,80,'b','filled');

s=trisurf(faces', vertices(1,:), vertices(2,:), vertices(3,:));
s.EdgeColor = [0,0,0];
s.FaceColor = 'none';
s=trisurf(faces_0', vertices_0(1,:), vertices_0(2,:), vertices_0(3,:)+1);
s.EdgeColor = [1,0,0];
s.FaceColor = 'none';
utils().apply_format(h,true,5);

if save_results
    utils().my_saveas(h,['graphics/grafica5_grid_',num2str(N),'iter_dots']);
    save(['graphics/grafica5_grid_',num2str(N),'iter'],'vertices','faces')
end