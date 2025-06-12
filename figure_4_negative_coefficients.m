% generates tikz figure

clear; close all; clc;

N = 7;
x = sym([0,-4,-3,-2,-1,4,3]');
y = sym([0,ones(1,4),1,-1]');
A = [ones(N,1) x y];

coef = zeros(1,N,'like',x);
for ind = 1:N
z = zeros(N,1);
z(ind) = 1;
cv = (A'*A)\(A'*z);
coef(ind) = [1,x(1),y(1)]*cv;
end

vertices = double([x,y,zeros(N,1)])';
faces = [
    1   1   1   1   1   1
    2   3   4   5   6   7
    3   4   5   6   7   2
];

% convert symbolic numbers to latex, including $ signs
labels = arrayfun(@(x) ['$' latex(x) '$'], coef, 'UniformOutput', false);
utils().make_draw('graphics/counterexample',vertices,faces,labels);

% subplot(1,2,1);
% plot(x,y,'ob');
% hold on;
% plot(x(ind),y(ind),'xr');

% subplot(1,2,2);
% plot3(x,y,z,'ob');
% hold on;
% plot3(x(ind),y(ind),z(ind),'xr');
% plot3(x,y,A*cv,'og');
% [X,Y] = meshgrid(linspace(-2,2,2));
% mesh(X,Y,cv(1)+cv(2)*X+cv(3)*Y)