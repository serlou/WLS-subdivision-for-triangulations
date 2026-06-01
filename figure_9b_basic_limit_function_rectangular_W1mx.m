% generates matlab figure

clear; close all;
load('datasets/large_rectangular.mat')
filename = 'graphics/grafica6_basic_rectangular_W1mx';
save_results = true;
L = 1.6;
weight = @(x) 1-x;

utils().basic_limit_function_figures(vertices, faces, L, weight, filename, save_results);