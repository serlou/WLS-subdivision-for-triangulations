function u = utils()
    % COMMON_PLOT_UTILS Provides common functions for visualizations
    %
    % This function returns a structure with common functions used
    % in different visualization techniques to avoid code duplication.
    %
    % Returns:
    %   utils: Structure with the following functions:
    %     - apply_format: Applies a common format to figures.
    %     - find_extraordinary_vertices: Finds vertices with valence different from 6.
    %     - generate_figures: Generates figures for the initial and final states of a mesh.
    %     - my_saveas: Saves figures in EPS format.
    %     - make_draw_ring: Creates a TikZ figure of a ring with vertices inside and outside.
    %     - make_figure_ring: Creates a MATLAB figure of a ring with vertices inside and outside.
    %     - new_figure: Creates a new figure with a specified position.
    %     - plotMesh: Plots a mesh given vertices and faces.
    %     - make_draw: Creates a TikZ figure of a mesh with labeled vertices.
    %     - basic_limit_function_figures: Generates figures for the basic limit function.
    
    
    u.apply_format = @apply_format;
    u.find_extraordinary_vertices = @find_extraordinary_vertices;
    u.generate_figures = @generate_figures;
    u.my_saveas = @my_saveas;
    u.make_draw_ring = @make_draw_ring;
    u.make_figure_ring = @make_figure_ring;
    u.new_figure = @new_figure;
    u.plotMesh = @plotMesh;
    u.make_draw = @make_draw;
    u.basic_limit_function_figures = @basic_limit_function_figures;
end

function apply_format(h,grid_on,num_grafica)
    if nargin < 2
        num_grafica = 0;
    end
    
    figure(h);
    colormap default
    axis tight;
    axis square;
    axis off;
    view(3);
    
    if grid_on
        h.Children(1).Children(1).LineWidth = 2;
    else
        h.Children(1).Children(1).LineStyle = 'none';
    end
    switch num_grafica
        case 1
        ll = light;
        ll.Position = [-1 -1 1];
        case 2
        ll = light;
        ll.Position = [-1 -1 1];
        axis([0.2 0.8 0.2 0.8])
        case 4
        ll = light;
        ll.Position = [-1 -1 1];
        axis([-2 2 -2 2])
        case 5
        axis([-1 1 -1 1])
        view(0,90)
        case 6
        ll = light;
        ll.Position = [-1 -1 1];
        case 7
        view(0,90)
        axis on;
        set(gca,'FontSize',20);
        % get axis limits
        xlim = get(gca,'XLim');
        ylim = get(gca,'YLim');
        % ceil them
        xlim = [floor(xlim(1))-1e-4, ceil(xlim(2))+1e-4];
        ylim = [floor(ylim(1))-1e-4, ceil(ylim(2))+1e-4];
        % set them back
        set(gca,'XLim',xlim);
        set(gca,'YLim',ylim);
        axis equal;
    end
end

function extraordinary_vertices = find_extraordinary_vertices(faces)
    % Find the list of vertices with valence different from 6
    % Input:
    %   faces - 3xN matrix with the faces of a triangulation
    % Output:
    %   extraordinary_vertices - list of vertices with valence different from 6
    
    % Find the number of vertices
    num_vertices = max(faces(:));
    
    % Initialize valence array
    valence = zeros(1, num_vertices);
    
    % Calculate the valence of each vertex
    for i = 1:size(faces, 2)
        for j = 1:3
            valence(faces(j, i)) = valence(faces(j, i)) + 1;
        end
    end
    
    % Find vertices with valence different from 6
    extraordinary_vertices = find(valence ~= 6);
end

function generate_figures(vertices, faces,vertices_0,faces_0,gn,save_results,filename,N)
    
    h = new_figure;
    scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),50,'filled');
    apply_format(h,true,gn);
    if save_results
        my_saveas(h,[filename,'_',num2str(0),'iter_dots']); %#ok<*UNRCH> 
    end
    
    h = new_figure;
    plotMesh(vertices_0, faces_0);
    apply_format(h,true,gn);
    if save_results
        my_saveas(h,[filename,'_',num2str(0),'iter']); %#ok<*UNRCH> 
    end
    
    h = new_figure;
    plotMesh(vertices, faces);
    apply_format(h,false, gn);
    if save_results
        my_saveas(h,[filename,'_',num2str(N),'iter']);
    end
    hold on;
    scatter3(vertices_0(1,:),vertices_0(2,:),vertices_0(3,:),50,'filled');
    if save_results
        my_saveas(h,[filename,'_',num2str(N),'iter_dots']);
    end
end

function my_saveas(h,nombre,pos)
    if nargin<3
        pos = [1           1        1440         900];
    end
    h.Position=pos;
    figure(h);
    drawnow;
    saveas(h,nombre,'epsc');
    % exportgraphics(h,[nombre,'.eps']);
end

function make_draw_ring(filename,vertices,faces,center, radius, eje)    
    if nargin < 6
        mini = min(vertices(1:2,:),[],2);
        width = max(vertices(2,:)) - min(vertices(2,:));
        width = max(1e-8,width);
    else
        mini = eje([1,3])';
        width = eje(4)-eje(3);
    end
    
    % get the vertex indices that are inside the ring
    distances = sqrt(sum((vertices(1:2,:) - center(1:2)).^2, 1));
    inside = distances < max(radius)-1e-15;
    
    fileID = fopen([filename,'.txt'],'w');
    fprintf(fileID,'\\begin{tikzpicture}\n');
    
    vertices = vertices(1:2,:);
    vertices = (vertices-mini)/width;
    vertex_in = vertices(:,inside);
    vertex_out = vertices(:,~inside);
    center = center(1:2);
    center = (center-mini)/width;
    radius = radius / width;

    % draw a green circle.
    fprintf(fileID,'\\draw[green] (%f,%f) circle (%f);\n',center(1),center(2),radius(1));
    if length(radius) > 1
        % draw a blue circle.
        fprintf(fileID,'\\draw[blue] (%f,%f) circle (%f);\n',center(1),center(2),radius(2));
    end
    
    order = [1:size(faces,1),1];
    edges = '\\draw[line width=0.05mm,   black] (%f,%f) -- (%f,%f);\n';
    for i=1:size(faces,2)
        for j=1:size(faces,1)
            fprintf(fileID,edges,...
            vertices(1,faces(order(j),i)),...
            vertices(2,faces(order(j),i)),...
            vertices(1,faces(order(j+1),i)),...
            vertices(2,faces(order(j+1),i)));
        end
    end
    
    points_red = '\\filldraw[red] (%f,%f) circle (0.1pt);\n';
    for i=1:size(vertex_in,2)
        fprintf(fileID,points_red,...
        vertex_in(1,i),vertex_in(2,i));
    end
    
    points_black = '\\filldraw[black] (%f,%f) circle (0.1pt);\n';
    for i=1:size(vertex_out,2)
        fprintf(fileID,points_black,...
        vertex_out(1,i),vertex_out(2,i));
    end
    
    % draw a red cross at the center
    fprintf(fileID,'\\draw[red,line width=0.1mm] (%f,%f) -- (%f,%f);\n',center(1)-0.02,center(2)-0.02,center(1)+0.02,center(2)+0.02);
    fprintf(fileID,'\\draw[red,line width=0.1mm] (%f,%f) -- (%f,%f);\n',center(1)-0.02,center(2)+0.02,center(1)+0.02,center(2)-0.02);
    
    fprintf(fileID,'\\end{tikzpicture}');
    fclose(fileID);
end

function make_figure_ring(vertices, faces, center, radius)
    figure;
    hold on;
    plotMesh(vertices,faces,3);
    view([0,90]);
    axis equal;
    axis off;
    
    % get the vertex indices that are inside the ring
    distances = sqrt(sum((vertices(1:2,:) - center(1:2)).^2, 1));
    inside = distances < max(radius)-1e-15;
    vertex_in = vertices(:,inside);
    vertex_out = vertices(:,~inside);
    
    
    % show dots at the vertices
    plot(vertex_in(1,:),vertex_in(2,:),'r.','MarkerSize',20);
    plot(vertex_out(1,:),vertex_out(2,:),'k.','MarkerSize',20);
    
    % mark the center with a red cross
    plot(center(1), center(2), 'rx', 'MarkerSize', 15, 'LineWidth', 5);
    
    % show green circle
    theta = linspace(0, 2*pi, 100);
    x = center(1) + radius(1) * cos(theta);
    y = center(2) + radius(1) * sin(theta);
    plot(x, y, 'g-', 'LineWidth', 2);
    
    if length(radius) > 1
        % show blue circle
        theta = linspace(0, 2*pi, 100);
        x = center(1) + radius(2) * cos(theta);
        y = center(2) + radius(2) * sin(theta);
        plot(x, y, 'b-', 'LineWidth', 2);
    end
end

function h = new_figure
    h = figure;
    % asd=get(groot,'MonitorPositions');
    % asd = asd(end,:);
    h.Position=[1           1        1440         900];
end

function plotMesh(vertices, faces,op)
    vertices=double(vertices);
    hold on;
    if nargin < 3
        op = 1; % default operation
    end
    switch op
        case 1
        s=trisurf(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        %             light('Position',[1 1 1])
        %             s.EdgeColor = 'none';
        s.EdgeColor = [0,0,0];
        s.FaceColor = 'interp';
        s.LineWidth = 2;
        %             colormap gray(1);
        colormap default
        case 2
        trimesh(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        colormap default
        %     colormap gray(1);
        case 3 % transparent faces
        s = trisurf(faces', vertices(1,:), vertices(2,:), vertices(3,:));
        s.EdgeColor = [0,0,0];
        s.LineWidth = 2;
        s.FaceColor = 'none';
    end
    axis tight;
    axis square;
    axis off;
    %     axis equal;
    view(3);
end

function make_draw(filename,vertices,faces,labels)
    % vertices contains the 2D coordinates
    % faces contains the indices of the vertices that are connected
    % labels is a cell containing string labels for each vertex
    mini = min(vertices,[],"all");
    maxi = max(vertices,[],"all");
    maxi = max(mini+1e-8,maxi);
    vertices = 10*(vertices-mini)/(maxi-mini);
    edges = '\\draw[line width=0.3mm,   black] (%.2f,%.2f) -- (%.2f,%.2f);\n';
    points = '\\filldraw[black] (%.2f,%.2f) circle (1.5pt);\n';
    % point labels where no point or edge is drawn
    points_labels = '\\node[anchor=south west] at (%.2f,%.2f) {%s};\n';
    
    fileID = fopen([filename,'.txt'],'w');
    fprintf(fileID,'\\begin{tikzpicture}\n');
    order = [1:size(faces,1),1];
    for i=1:size(faces,2)
        for j=1:size(faces,1)
            fprintf(fileID,edges,...
            vertices(1,faces(order(j),i)),...
            vertices(2,faces(order(j),i)),...
            vertices(1,faces(order(j+1),i)),...
            vertices(2,faces(order(j+1),i)));
        end
    end
    for i=1:size(vertices,2)
        fprintf(fileID,points,...
        vertices(1,i),vertices(2,i));
        fprintf(fileID,points_labels,...
        vertices(1,i),vertices(2,i),labels{i});
    end
    fprintf(fileID,'\\end{tikzpicture}');
    fclose(fileID);
end

function basic_limit_function_figures(vertices, faces, L, weight, filename, save_results)
    % find the [0;0;0] vertex
    [~, idx] = find(vertices(1,:)==0 & vertices(2,:)==0 & vertices(3,:)==0);
    vertices(3,idx) = 1;
    
    vertices_0 = vertices;
    faces_0 = faces;
    
    N=5;
    for i=1:N
        [vertices, faces] = RegressionSubdivision(vertices, faces, L/2^(i-1), weight,'functional');
    end 
    
    
    h = utils().new_figure();
    utils().plotMesh(vertices_0, faces_0);
    utils().apply_format(h,true,6);
    if save_results
        utils().my_saveas(h,[filename,'_',num2str(0),'iter']); %#ok<*UNRCH> 
    end

    bool = vertices(3,:) > 0;
    vertices(3,~bool) = NaN;
    
    h = utils().new_figure();
    utils().plotMesh(vertices, faces);
    utils().apply_format(h,false,6);
    if save_results
        utils().my_saveas(h,[filename,'_',num2str(N),'iter']);
    end
    save([filename,'_',num2str(N),'iter'],'vertices','faces')
    
    
    h = utils().new_figure();
    utils().plotMesh(vertices, faces);
    utils().apply_format(h,false,7);
    if save_results
        utils().my_saveas(h,[filename,'_support']);
    end
end