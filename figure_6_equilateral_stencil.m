% generates matlab figure
% generates tikz figure
clear;
close all;

load('datasets/equilateral_for_figure.mat');
radius = [3/2,sqrt(3)];

centers = [0, 0, 0;
1/2, 0, 0;
1/4, sqrt(3)/4, 0;
3/4, sqrt(3)/4, 0]';

ori = 'rhvd';

labels{1}={'c','c','c','f','c','c','c'};
labels{2}={'b','d','b','a','e','e','a','b','d','b'};
labels{3}={'a','b','b','e','d','d','e','b','b','a'};
labels{4}={'b','a','d','e','b','b','e','d','a','b'};

for i=1:length(labels)
    for j=1:length(labels{i})
        labels{i}{j} = ['$',labels{i}{j},'$'];
    end
end

for i = 1:size(centers, 2)
    center = centers(:, i);
    utils().make_figure_ring(vertices, faces, center, radius,labels{i});
    utils().make_draw_ring(['graphics/equilateral_ring_',ori(i)], vertices, faces, center, radius,labels{i});
end