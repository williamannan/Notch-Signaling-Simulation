function [Ang] = FilopAngles(Nc, r, FixedPoint)

% obtain the vertex and the center coordinate for each cell

[VertCod,CenterCod, ~] = TwoDGeom(Nc,r);

% FixedPoint = ['NE', 'NW', 'NW', 'SE', 'C']
% Remember to specify only one of the specified fixed location

if strcmp(FixedPoint,'NE')
    cellNum = Nc*Nc;
    P = VertCod(cellNum, :, 1);

elseif strcmp(FixedPoint,'NW')
    cellNum = Nc*(Nc - 1) + 1;
    P = VertCod(cellNum, :, 3);

elseif strcmp(FixedPoint,'SW')
    cellNum = 1;
    P = VertCod(cellNum, :, 4);

elseif strcmp(FixedPoint,'SE')
    cellNum = Nc;
    P = VertCod(cellNum, :, 6);

elseif strcmp(FixedPoint,'C')
    if mod(Nc, 2)==0
        cellNum = (Nc/2)*(Nc-1);
        P = VertCod(cellNum, :, 2);
    else
        cellNum = Nc*((Nc-1)/2) + (Nc+1)/2;
        P = CenterCod(cellNum, :);
    end
end

%%%%----------------

Ang = zeros(Nc*Nc, 6);  %% initialize angle matrix

for icell = 1: Nc*Nc
    % Extracting the vertices for any given cell
    V = squeeze(VertCod(icell, :,:))';   
    
    % Compute the corresponding edge vector
    d = [V(6,:)-V(1,:);     %% Creating the edge vectors
        V(1,:)-V(2,:);
        V(2,:)-V(3,:);
        V(3,:)-V(4,:);
        V(4,:)-V(5,:);
        V(5,:)-V(6,:)];
    % Compute the vector connecting the vertex to the fixed point
    VP = P - V;

    % Compute the angles VP and d make with the positive x-axis
    theta_d = atan2(d(:,2),d(:,1));
    theta_vp = atan2(VP(:,2),VP(:,1));
    
    % Transform the angle to clockwise rotation and store
    Ang(icell, :) = mod(theta_vp-theta_d, 2*pi);

end


end


