function A = compute_mask(w,L,L_stencil,basis,draw)
if nargin < 5
    draw = false;
end

ind = [
    0,0
    1,0
    0,1
    1,1]'; % for all uniform grids

insertion_point = ind.*0.5; % for all uniform grids

of = ceil(4*L_stencil);
A = zeros(2*of,2*of,'like',basis);

min1 = 2*of+1;
min2 = 2*of+1;
max1 = 0;
max2 = 0;

close all;
for j = 1:length(insertion_point)

    [stencil,a] = define_weighted(insertion_point(:,j),w,L,L_stencil,basis,draw);

    for i = 1:length(stencil)
        i1 = of-2*stencil(1,i)+ind(1,j);
        i2 = of-2*stencil(2,i)+ind(2,j);
        A(i1,i2) = a(i);
        min1 = min(min1,i1);
        min2 = min(min2,i2);
        max1 = max(max1,i1);
        max2 = max(max2,i2);
    end
end
% remove those rows and columns that are all zeros

if isa(A, 'sym')
    A = simplify(A(min1:max1,min2:max2));
end
end

function [stencil_0,a] = define_weighted(insertion_point,w,L,L_stencil,base,draw)
if ~iscolumn(insertion_point)
    insertion_point = insertion_point';
end
of = ceil(8*L_stencil);
[X,Y]=meshgrid(-of:of,-of:of);
stencil_0 = transpose([X(:),Y(:)]);
base = transpose(base);
stencil = stencil_0 - insertion_point;
stencil = base*stencil;
insertion_point = [0,0]';
inside_ball = sum(stencil.^2,1) < L_stencil^2;
stencil = stencil(:,inside_ball);
stencil_0 = stencil_0(:,inside_ball);

n = size(stencil,2);

W = zeros(n,n,'like',base);
for i=1:n
    W(i,i) = w(norm(stencil(:,i))/L);
end

if draw
    figure('Position',[0 0 1200 800]);
    plot(stencil(1,:),stencil(2,:),'*')
    hold on;
    plot(insertion_point(1),insertion_point(2),'ok','MarkerSize',10)
end

% Method 1 (Fast): Since the grid is uniform, using formula \alpha_i = w_i / sum_j w_j
a = diag(W);
a = transpose(a) / sum(a);

% Method 2: Solving a linear system for each coefficient
% a_comparison = a;
% a = sym(zeros(1,n));
% for i = 1:n
%     data = sym(zeros(1,n));
%     data(i) = 1;

%     A = sym([ones(size(stencil,2),1),stencil']);

%     coef = transpose(A)*W*A\(transpose(A)*W*data');
%     a(i) = simplify(coef(1));
%     if draw
%         if all(isreal(stencil))
%             text(stencil(1,i),stencil(2,i),['$$',latex(a(i)),'$$'],'Interpreter', 'latex', 'FontSize', 8);
%         end
%     end
% end
% % assert that both methods give the same result
% assert(isequal(simplify(a - a_comparison),zeros(1,n,'like',a)));

if draw
    axis(double([-L_stencil-1,L_stencil+1,-L_stencil-1,L_stencil+1]))
    axis equal;
end
end