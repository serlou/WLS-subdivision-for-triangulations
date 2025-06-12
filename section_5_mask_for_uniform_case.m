close all;

% w = @(x) 1;
w = @(x) 1 - x;
% syms W(x)
% w = @(x) W(x);

syms L  % actual value for the mask coefficient

% 3/2 < 1.6 < sqrt(3) < sqrt(13)/2
L_stencil = 1.6; % fake value just to compute the stencil

% rectangular grid:
basis = sym([1,0; 0,1]);

% equilateral grid:
% basis = [1,0; 0.5,sqrt(sym(3))/2];

A = compute_mask(w,L,L_stencil,basis);

pretty(A);
latex(A)

AA = simplify(72*A);
pretty(AA);
latex(AA)