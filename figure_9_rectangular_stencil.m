% generates matlab figure
% generates tikz figure
clear;
close all;

load('datasets/rectangular_for_figure.mat');
radius = [3/2,sqrt(13)/2];

centers = [0, 0, 0;
1/2, 0, 0;
0, 1/2, 0;
1/2, 1/2, 0]';

ori = 'rhvd';

labels{1}={'c','e','c','e','h','e','c','e','c'};
labels{2}={'d','d','b','g','g','b','d','d'};
labels{3}={'b','d','g','d','d','g','d','b'};
labels{4}={'a','a','a','f','f','a','a','f','f','a','a','a'};

for i=1:length(labels)
    for j=1:length(labels{i})
        labels{i}{j} = ['$',labels{i}{j},'$'];
    end
end

for i = 1:size(centers, 2)
    center = centers(:, i);
    utils().make_figure_ring(vertices, faces, center, radius, labels{i});
    utils().make_draw_ring(['graphics/rectangular_ring_',ori(i)], vertices, faces, center, radius, labels{i});
end