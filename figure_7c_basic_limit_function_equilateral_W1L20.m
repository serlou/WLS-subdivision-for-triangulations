clear;
close all;
save_results = true;
filename = 'graphics/grafica6_basic_equilateral_W1_L20';

% equilateral grid:
basis = [1,0; 0.5,sqrt(3)/2];

L = 20;
W = @(x) 1;

utils().basic_limit_function_figures_fast(W, L, basis, filename, save_results);