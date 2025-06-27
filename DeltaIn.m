function [DinJ, DinF, FCount] = DeltaIn(Nc, r, F, Ang, CellsAway,a_apical, a_basal, Delta)
[DirVec, Fbase,FilopDisc, JunctDisc,~] = Filop_vectors(Nc, r, F, Ang, CellsAway);

DinJ = zeros(1, Nc*Nc);
DinF = zeros(1, Nc*Nc);
FCount = zeros(1, Nc*Nc);

for cellp = 1:Nc*Nc
    SumOfDeltaJ = 0;
    SumOfDeltaF = 0;
    count = 0;
    %%% -----Conputing Din for junctional interaction---
    JList = JunctDisc{cellp};
    for ijunc = 1:length(JList)
        SumOfDeltaJ = SumOfDeltaJ + (a_apical/2)*Delta(JList(ijunc));
    end
    DinJ(cellp) = SumOfDeltaJ;


    %%% -----Conputing Din for junctional interaction--------
    FList = FilopDisc{cellp};
    for ifilop = 1:length(FList)
        for ip = 1:6
            Vp = DirVec(cellp, :, ip);
            bp = Fbase(cellp, :, ip);

            for jq = 1:6
                Vq = DirVec(FList(ifilop), :, jq);
                bq = Fbase(FList(ifilop), :, jq);

                A = [Vp', -Vq'];
                bb = bq' - bp';
                if abs(det(A)) > 1e-10  %% Making sure the matrix A is non-singular
                    T = A\bb;          %% Computing T
                    if ( T(1)>=0 && T(1) <= 1 ) && ( T(2)>=0 && T(2) <= 1 )
                        SumOfDeltaF = SumOfDeltaF + (a_basal/24)*(Delta(FList(ifilop)));
                        count = count +1;
                        break
                    end

                end
            end

        end

    end

    DinF(cellp) = SumOfDeltaF;
    FCount(cellp) = count;
end

%%%%%%%%%%%%%%%%%%%%%%========================================================================
%% -------- Section below is just for visualization, comment it out during simulation --------
%%%%%%%%%%%%%%%%%%%%%%========================================================================
figure
DinTot = DinJ + DinF;
cellnum = 1:Nc*Nc;
yyaxis left
plot(cellnum,DinTot, 'c -', cellnum,DinJ, 'm.-', cellnum,DinF, 'b :', 'LineWidth', 1.5);
title('Din profile')
xlabel('Cell number', 'FontSize',13);
ylabel('Din', 'FontSize',13);
ax = gca;                    % Get current axes
ax.YAxis(1).Color = 'k';     % Set left y-axis color to blue

yyaxis right
plot(cellnum, FCount, 'r.-', 'LineWidth',1.5)
ylabel('Filopodia count', 'FontSize',13);

ax = gca;                    % Get current axes
ax.YAxis(2).Color = 'r';     % Set left y-axis color to blue

legend('Total-Din', 'Junctional Din', 'Filopodia Din', '# Filopodia Intersection','Location','north');

%%% ------------- Making a heatmap for the cells --------------
[~,~, Pcell] = TwoDGeom(Nc,r);
figure;
hold on;
axis equal;
colormap(jet);  % You can change this to parula, hot, etc.

% Normalize Fcount to map to color indices
Fmin = min(FCount);
Fmax = max(FCount);
C = (FCount - Fmin) / (Fmax - Fmin);    % Scale to [0,1]
colors = colormap;                      % Get current colormap
nColors = size(colors, 1);

for i = 1:length(Pcell)
    if ~isempty(Pcell(i).Vertices)
        % Find color index based on normalized value
        idx = max(1, round(C(i) * (nColors - 1)) + 1);
        fill(Pcell(i).Vertices(:,1), Pcell(i).Vertices(:,2), colors(idx, :), 'EdgeColor', 'k');
    end
end

colorbar;
clim([Fmin Fmax]);                           % Make colorbar match Fcount range
title('Heatmap of Polygon Intersections');

end