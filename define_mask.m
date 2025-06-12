close all;
A = sym([]);

w = @(x) 1;
% w = @(x) 1 - x;
% syms phi(x)
% w = @(x) phi(x);

% syms L
L = 2;

% general grid:
% syms v1 [2,1]
% syms v2 [2,1]

% rectangular grid:
% base = sym([1,0; 0,1]);

% equilateral grid:
base = [1 1;0 1]/[1 cos(sym(pi)/3);0 sin(sym(pi)/3)];

insertion_point = [
    0,0
    1,0
    0,1
    1,1]'.*0.5; % for all uniform grids

ind = [
    -1 -1
    -1 0
    0 -1
    0 0];

close all;
for j = 1:length(insertion_point)
    if islogical(base) || all(isreal(base))
        figure('Position',[0 0 1200 800])
    end
    [stencil,a] = utils().define_weighted(insertion_point(:,j),w,L,base);
    stencil = stencil + (1+ind(j,:)')/2 + 2;
    for i = 1:length(stencil)
        A(2*stencil(1,i)-ind(j,1),2*stencil(2,i)-ind(j,2)) = a(i);
    end
    axis([-2,3,-2,3])
    axis equal
end
flipud(simplify(A(2:end,2:end)))
lat = latex(ans);
lat2 = strrep(lat,'v_{11}','v^1_1');
lat2 = strrep(lat2,'v_{21}','v^2_1');
lat2 = strrep(lat2,'v_{12}','v^1_2');
lat2 = strrep(lat2,'v_{22}','v^2_2')
pretty(flipud(simplify(A(2:end,2:end))))