% generates matlab figure

clear; close all;
load('datasets/dataset2.mat')
save_results = true;
weight = @(x) 1-x;
L = 1;

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
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',num2str(0),'iter']); %#ok<*UNRCH> 
end
h = utils().new_figure();
utils().plotMesh(vertices, faces);
utils().apply_format(h,false,4);
if save_results
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',num2str(N),'iter']);
end
hold on;
scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),25,'filled');
if save_results
    utils().my_saveas(h,['graphics/grafica4_sin_cos_',num2str(N),'iter_dots']);
    save(['graphics/grafica4_sin_cos_',num2str(N),'iter'],'vertices','faces')
end

format short e
x_int = [-1,1];
y_int = [-1,1];
bool_int = (vertices(1,:)>=x_int(1)) & (vertices(1,:)<=x_int(2)) & (vertices(2,:)>=y_int(1)) & (vertices(2,:)<=y_int(2));
MAE_int = norm(vertices(3,bool_int)-F(vertices(1,bool_int),vertices(2,bool_int)),'Inf')/length(vertices(3,bool_int))
MSE_int = norm(vertices(3,bool_int)-F(vertices(1,bool_int),vertices(2,bool_int)),2)/length(vertices(3,bool_int))

initial_MAE = norm(vertices_0(3,:)-F(vertices_0(1,:),vertices_0(2,:)),'Inf')/length(vertices_0(3,:))
initial_MSE = norm(vertices_0(3,:)-F(vertices_0(1,:),vertices_0(2,:)),2)/length(vertices_0(3,:))

format default