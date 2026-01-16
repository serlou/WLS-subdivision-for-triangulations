% TODO: Dioni, escribe aqui el código necesario para generar la nueva tabla:
clear;

L = 5;

% rectangular grid:
basis = [1,0; 0,1];

% equilateral grid:
% basis = [1,0; 0.5,sqrt(sym(3))/2];

v(1)=double(utils().get_noise_reduction_coefficient(@(x) 1,L,L,basis));
v(2)=double(utils().get_noise_reduction_coefficient(@(x) 1-x,L,L,basis));
v(3)=double(utils().get_noise_reduction_coefficient(@(x) 1-x.^2,L,L,basis));
v(4)=double(utils().get_noise_reduction_coefficient(@(x) (1-x.^2).^2,L,L,basis));
v(5)=double(utils().get_noise_reduction_coefficient(@(x) (1-x.^3).^3,L,L,basis));
v(6)=double(utils().get_noise_reduction_coefficient(@(x) (1-x.^2).^3,L,L,basis));
v(7)=double(utils().get_noise_reduction_coefficient(@(x) exp(-x),L,L,basis));

v'