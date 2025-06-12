% generates matlab figure
% generates tikz figure

load('datasets/triangulation2.mat');
center = vertices(:,6);
radius = 1.4;
utils().make_figure_ring(vertices, faces, center, radius);
utils().make_draw_ring('graphics/ring_replacement', vertices, faces, center, radius);