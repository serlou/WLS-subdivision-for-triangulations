function [V, F] = RegressionSubdivision (vertices, faces, L, weight, mode)
    % RegressionSubdivision Subdivides a mesh using regression-based rules.
    % Warning: The implementation is far from optimal and is intended
    % for research purposes. For large amount of vertices, it is
    % memory inefficient.
    % 
    % Inputs:
    %   vertices - 3xN matrix of vertex coordinates
    %   faces - 3xM matrix of face indices
    %   L - scalar parameter controlling the subdivision stencil size
    %   mode - string specifying the mode ('functional' or 'svd')
    %
    % Outputs:
    %   V - updated vertices after subdivision
    %   F - updated faces after subdivision
    
    vertices_update = zeros(size(vertices));
    new_faces = zeros(3, size(faces, 2) * 4);
    countNewF = 0;
    new_vertices = zeros(3, size(faces, 2) * 3);
    new_VMatrix = sparse(size(vertices, 2), size(vertices, 2));
    count_newV = 0;
    count_initialV = size(vertices,2);
    v_neighs = preprocessing(faces, count_initialV);
    
    if nargin < 4
        mode = 'functional';
    end
    switch mode
        case 'svd'
        rule = @(v,center) RuleSVD(v,center,L,weight);
        coord = 1:3;
        otherwise
        rule = @(v,center) Rule(v,center,L,weight);
        coord = 1:2;
    end
    
    order = [
    1 2 3
    3 1 2
    2 3 1];
    
    for f = 1 : size(faces,2)
        v_new = zeros(1,3);
        for i= 1:3
            v(1) = faces(order(i,1), f);
            v(2) = faces(order(i,2), f);
            v(3) = faces(order(i,3), f);
            
            if new_VMatrix(v(1), v(2)) == 0
                % If the distance between vertices v1 and v2 is greater than or equal to 2*L, throw an error
                if norm(vertices(coord,v(1))-vertices(coord,v(2))) >= 2*L
                    % With this L value, the stencil is too small to include the vertices v1 and v2
                    figure(99); clf;
                    plot3(vertices(1,:),vertices(2,:),vertices(3,:),'.');
                    hold on;
                    plot3(vertices(1,v),vertices(2,v),vertices(3,v),'r-');
                    error('L must be greater');
                end
                
                %% insertion
                mid_point = mean(vertices(:,v(1:2)),2);
                stencil = compute_ring(v(1:2),mid_point(coord),L,vertices(coord,:),v_neighs);
                q = rule(vertices(:,stencil),mid_point);
                count_newV = count_newV + 1;
                new_vertices( :, count_newV) = q;
                new_VMatrix(v(1), v(2)) = count_newV;
                new_VMatrix(v(2), v(1)) = count_newV;
                v_new(i) = count_newV;
            else
                v_new(i) = new_VMatrix(v(1), v(2));
            end
        end
        
        countNewF = countNewF + 1;
        new_faces(1, countNewF) = v(1);
        new_faces(2, countNewF) = v_new(1) + count_initialV;
        new_faces(3, countNewF) = v_new(3) + count_initialV;
        
        countNewF = countNewF + 1;
        new_faces(1, countNewF) = v(2);
        new_faces(2, countNewF) = v_new(2) + count_initialV;
        new_faces(3, countNewF) = v_new(3) + count_initialV;
        
        
        countNewF = countNewF + 1;
        new_faces(1, countNewF) = v(3);
        new_faces(2, countNewF) = v_new(1) + count_initialV;
        new_faces(3, countNewF) = v_new(2) + count_initialV;
        
        countNewF = countNewF + 1;
        new_faces(1, countNewF) = v_new(1) + count_initialV;
        new_faces(2, countNewF) = v_new(2) + count_initialV;
        new_faces(3, countNewF) = v_new(3) + count_initialV;
    end
    
    for i = 1 : size(vertices,2)
        stencil = compute_ring(i,vertices(coord,i),L,vertices(coord,:),v_neighs);
        % disp(['Cantidad de puntos en el stencil: ',num2str(length(stencil))]);
        vertices_update(:,i) = rule(vertices(:,stencil),vertices(:,i));
    end
    
    new_vertices = new_vertices(:,1:count_newV);
    V = [vertices_update new_vertices];
    F = new_faces;
end

function [q] = RuleSVD(v,center,L,weight)
    v = v-center;
    [U,~,~] = svd(v);
    v_svd = U\v;
    % plot3(v_svd(1,:),v_svd(2,:),v_svd(3,:),'.');
    % assert(all(abs(sqrt(sum(v_svd.^2,1))-sqrt(sum(v.^2,1)))<1e-14));    
    
    % dist = sqrt(sum((v(1:2,:)).^2,1));
    % assert(all(dist<L));
    
    [q] = Rule(v_svd,[0;0;0],L,weight);
    q = U*q + center;
end

function [q] = Rule(v,center,L,weight)
    switch size(v,2)
        case 0
        error('L is too small');
        case 1
        if norm(v(1:2)-center(1:2),'Inf') < 1e-14
            q = center;
            return;
        else
            error('L is too small');
        end
    end
    
    w = compute_weights(v,center,L,weight);
    q = zeros(3,1);
    q(1:2) = center(1:2);
    v(1:2,:) = v(1:2,:) - center(1:2);
    q(3) = perform_regression(v,w);
end

function [qx] = perform_regression(v,w)
    W = sum(w);
    switch size(v,2)
        case 1
        qx = v(3);
        return;
        case 2
        qx = v(3,:)*w'/W;
        return;
    end
    
    X = sum(w.*v(1,:));
    X2 = sum(w.*v(1,:).^2);
    Y = sum(w.*v(2,:));
    Y2 = sum(w.*v(2,:).^2);
    XY = sum(w.*v(1,:).*v(2,:));
    F = sum(w.*v(3,:));
    FX = sum(w.*v(3,:).*v(1,:));
    FY = sum(w.*v(3,:).*v(2,:));
    A = [
    W   X   Y
    X   X2  XY
    Y   XY  Y2
    ];
    A1 = [
    F   X   Y
    FX   X2  XY
    FY   XY  Y2
    ];
    
    detA = det(A);
    if abs(detA)<1e-14
        error('Three points of a face are almost collinear. Points: [%f, %f, %f]', v(1), v(2), v(3));
    end
    
    qx = det(A1)/detA;
end

function w = compute_weights(v,center,L,weight)
    dist = sqrt(sum((v(1:2,:) - center(1:2)).^2,1));
    % assert(all(dist<L));
    w = weight(dist/L);
    w(dist/L>=1) = 0;
    w = w/sum(w);
end

function [v_neighs] = preprocessing(faces, n_vertices)
    % This function could be much faster if it is implemented in C++, using a
    % std::vector<std::set<int>> to store v_neighs.
    
    v_neighs = zeros(n_vertices,1);
    
    for f = 1: size(faces,2)
        v(1) = faces(1, f);
        v(2) = faces(2, f);
        v(3) = faces(3, f);
        
        v_neighs(v(1),1) = v_neighs(v(1),1) + 2;
        v_neighs(v(2),1) = v_neighs(v(2),1) + 2;
        v_neighs(v(3),1) = v_neighs(v(3),1) + 2;
        v_neighs(v(1),v_neighs(v(1),1)+(0:1)) = v([2,3]);
        v_neighs(v(2),v_neighs(v(2),1)+(0:1)) = v([1,3]);
        v_neighs(v(3),v_neighs(v(3),1)+(0:1)) = v([1,2]);
    end
    for i = 1:n_vertices
        neighs = unique(v_neighs(i,2:v_neighs(i,1)+1));
        v_neighs(i,1) = length(neighs);
        v_neighs(i,2:v_neighs(i,1)+1) = neighs;
        v_neighs(i,v_neighs(i,1)+2:end) = 0;
    end
    v_neighs = v_neighs(:,1:max(v_neighs(:,1))+1);
end

% function [v_neighs] = preprocessing2(faces, n_vertices)
% It is equivalent to the previous function, but it is faster
% when the number of vertices is small.
%     v_connections = zeros(n_vertices, n_vertices);

%     for f = 1: size(faces,2)
%         v(1) = faces(1, f);
%         v(2) = faces(2, f);
%         v(3) = faces(3, f);

%         v_connections(v(1), v(2)) = v(3);
%         v_connections(v(2), v(1)) = v(3);
%         v_connections(v(1), v(3)) = v(2);
%         v_connections(v(3), v(1)) = v(2);
%         v_connections(v(2), v(3)) = v(1);
%         v_connections(v(3), v(2)) = v(1);
%     end
%     v_neighs = zeros(n_vertices,0);
%     for i = 1:n_vertices
%         neighs = find(v_connections(i,:));
%         v_neighs(i,1) = length(neighs);
%         v_neighs(i,2:v_neighs(i,1)+1) = neighs;
%     end
% end

function neigh = find_neighbours(v1, v_neighs)
    neigh = v_neighs(v1,2:v_neighs(v1,1)+1);
end

function stencil = compute_ring(stencil, center, L, vertices, v_neighs, stencil_obtained)
    % This function could be much faster if it is implemented in C++, using a
    % std::set to store the stencil.
    
    if nargin < 6
        stencil_obtained = [];
    end
    to_be_added = [];
    for v1=stencil
        neigh = find_neighbours(v1, v_neighs);
        vertices_in_range = sqrt(sum((vertices(:,neigh)-center).^2,1))<L;
        to_be_added = [to_be_added, neigh(vertices_in_range)]; %#ok<AGROW>
    end
    stencil = [stencil_obtained,stencil];
    to_be_added = setdiff(to_be_added,stencil);
    if ~isempty(to_be_added)
        stencil = compute_ring(to_be_added, center, L, vertices, v_neighs, stencil);
    end
    % assert(all(sqrt(sum((vertices(:,stencil)-center).^2,1))<L));
end