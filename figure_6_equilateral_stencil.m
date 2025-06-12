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

for i = 1:size(centers, 2)
    center = centers(:, i);
    utils().make_figure_ring(vertices, faces, center, radius);
    utils().make_draw_ring(['graphics/equilateral_ring_',ori(i)], vertices, faces, center, radius);
end