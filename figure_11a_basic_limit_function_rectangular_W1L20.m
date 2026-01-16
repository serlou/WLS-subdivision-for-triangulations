clear;
close all;
save_results = true;
filename = 'graphics/grafica6_basic_rectangular_W1_L20';

% rectangular grid:
basis = [1,0; 0,1];

L = 20;
W = @(x) 1;

utils().basic_limit_function_figures_fast(W, L, basis, filename, save_results);