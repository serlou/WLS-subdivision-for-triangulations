% generates matlab figure

clear; close all;
load('datasets/sphere100.mat')
filename = 'graphics/grafica1_sphere';
save_results = true;

% random seed
rng(1);

vertices = vertices_0;
vertices = vertices + 0.05*(randn(1,size(vertices,2))).*vertices;
vertices_0 = vertices;
faces_0 = faces;
weight = @(x) 1-x;
L = 0.6;

N=5;
for i=1:N
    [vertices, faces] = RegressionSubdivision(vertices, faces, L, weight,'svd');
    L = L/2;
end
utils().generate_figures(vertices, faces,vertices_0,faces_0,1,save_results,filename,N);
save([filename,'_',num2str(N),'iter'],'vertices','faces')