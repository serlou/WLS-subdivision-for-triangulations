% TODO: Dioni, escribe aqui el código necesario para generar la nueva tabla:

L = 1.6;

% rectangular grid:
basis = sym([1,0; 0,1]);

% equilateral grid:
% basis = [1,0; 0.5,sqrt(sym(3))/2];

double(get_noise_reduction_coefficient(@(x) 1,L,L,basis))
double(get_noise_reduction_coefficient(@(x) 1-x,L,L,basis))