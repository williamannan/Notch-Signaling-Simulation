function [VertCod,CenterCod, Pcell] = TwoDGeom(Nc,r)

%% Creating the vertices of the first cell centered at the origin.
% theta = (pi/6)*[1,3,5,7,9,11];  %% Angle of lines connection vertices to the center
kk = 0:5;
theta = (pi/6)*(2*kk+1);
V(1,:) = r*cos(theta);           %% x-components of one cell
V(2,:) = r*sin(theta);           %% y-components of one cell

%% Creating the coordinate of the cell centers
xCent= zeros(Nc*Nc,1);           %% Initializing matrix for the x-component of cell centers
yCent= zeros(Nc*Nc,1);           %% Initializing matrix for the y-component of cell centers

%%% ----x-components of the centers for the first row of cells -------
row1 = zeros(Nc,1);              
for j = 2:Nc
    row1(j) = row1(j-1) + 2*r*cos(pi/6);
end
xCent(1:Nc) = row1;  

%%% ----x and y-components of the centers for the cells in the tissue -------
for irow = 2:Nc
    currentRow = (irow-1)*Nc+1:irow*Nc;        %% index for the row of cells under consideration
    PreviousRow = (irow-2)*Nc+1:(irow-1)*Nc;   %% index for the previous row of cells 
    
    %%% ----Computing the y coordinates for the cell centers --------
    yCent(currentRow) = yCent(PreviousRow) + r*(1+sin(pi/6));

    %%% ----Computing the x coordinates for the 2nd to last row of cell centers --------
    if mod(irow,2) == 0
        xCent(currentRow) = row1 - r*cos(pi/6);
    else
        xCent(currentRow) = row1;
    end

end

CenterCod = [xCent, yCent];             %% (x,y)-Coordinates at the cell centers

%% Creating the hexagonal polygons as well as the coordiate at the vertices
VertCod = zeros(Nc*Nc, 2,6);            %% Coordinate at the vertices 
Pcell = repmat(polyshape(), 1, Nc*Nc);  %% Predefining polyshape object

for jcell = 1:Nc*Nc
    Vert = V + CenterCod(jcell, :)';    %%translating the first hexagon to any location in the tissue
    VertCod(jcell, :,:) = Vert;         

    Pcell(jcell) =polyshape(Vert');     %% Polygonal shapes for the cells
end

%%%%%%%%%==============================================================
%%%---------- visualizing the tissue ----------------------
figure
plot(Pcell,'FaceColor','c','edgecolor','k','facealpha',1);
% title('\bf Visualizing the geometry', 'FontSize',13)
daspect([1, 1, 1])
axis tight
box on
pg = gcf;
exportgraphics(pg,'2D_Geom.png','Resolution',300);


% % set(gca, 'XTick', [], 'YTick', []);     % Remove ticks
% % set(gca, 'XColor', 'k', 'YColor', 'k'); % Optional: keep axis lines black
% % box on; 
% % axis tight
% % pg = gcf;
% % exportgraphics(pg,'2D_Geom.png','Resolution',300);

end
