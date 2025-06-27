function [DirVec, Fbase,FilopDisc, JunctDisc,Ftip] = Filop_vectors(Nc, r, F, Ang, CellsAway)

%%% [Vertex, Center, Pcell] = VisualTwoDGeomTwoDGeom(Nc,r);
[VertCod,CenterCod, Pcell] = TwoDGeom(Nc,r);


Fbase = VertCod;               %% Coordinate at the base of the filopodia (vertices of the cell)
Ftip = zeros(Nc*Nc,2,6);       %% Coordinate at the end of the filopodia
DirVec = zeros(Nc*Nc,2,6);     %% Direction vector along the filopodia

for jcell = 1:Nc*Nc            %% looping over the cells
    for jvert = 1:6
        if jvert == 1
            a_end = Fbase(jcell, :, 6);
            a_start = Fbase(jcell,:,jvert);
            vhat = (a_end - a_start)/(norm(a_end - a_start));  %% unit vector
            R = [cos(Ang(jcell,jvert)), -sin(Ang(jcell,jvert)); sin(Ang(jcell,jvert)), cos(Ang(jcell,jvert))];

            u = F(jcell, jvert)*R*vhat';           %% vector in the direction of the filopodia

            DirVec(jcell,:,jvert) = u;                %% vector in the direction of the filopodia

            Ftip(jcell, :, jvert) = a_start + u';
        else
            a_end = Fbase(jcell, :, jvert-1);
            a_start = Fbase(jcell,:,jvert);
            vhat = (a_end - a_start)/(norm(a_end - a_start));  %% unit vector
            R = [cos(Ang(jcell,jvert)), -sin(Ang(jcell,jvert)); sin(Ang(jcell,jvert)), cos(Ang(jcell,jvert))];

            u = F(jcell, jvert)*R*vhat';        %% vector in the direction of the filopodia

            DirVec(jcell,:,jvert) = u;                %% vector in the direction of the filopodia

            Ftip(jcell, :, jvert) = a_start + u';
        end
    end
end

%% Creating the contact list
FilopDisc = cell(Nc*Nc, 1);
JunctDisc = cell(Nc*Nc, 1);

for icell = 1: Nc*Nc
    Fcount = 1;
    Jcount = 1;
    FList = zeros(1, Nc*Nc); %[];
    JList = zeros(1, Nc*Nc); %[];

    for rcell = 1:Nc*Nc
        if ( icell ~= rcell)&&( norm(CenterCod(icell, :) - CenterCod(rcell, :)) <= 2*r )
            JList(Jcount) = rcell;
            Jcount = Jcount+1;
        end

        if (icell ~= rcell) && (norm( CenterCod(icell, :) - CenterCod(rcell, :) ) <= 2*r*CellsAway)
            FList(Fcount) = rcell;
            Fcount = Fcount+1;
        end
    end

    FilopDisc{icell} = FList(1:Fcount-1);
    JunctDisc{icell} = JList(1:Jcount-1);
end

%%%%%%%%=============================================
%%%%-----Visualize the output ----------------------
prob = 0.01;                                  % Desired proportion of 1s
kprob = round(prob*Nc*Nc);                   % Number of 1s
V = [ones(1, kprob), zeros(1, Nc*Nc - kprob)];  % Create vector with k ones and N-k zeros
Epit_SOP = V(randperm(Nc*Nc))';

figure
plot(Pcell,'FaceColor','c','edgecolor','k','facealpha',1);
hold on
for mm = 1:Nc*Nc
    % if Epit_SOP(mm) == 1
    %     plot(Pcell(mm),'FaceColor','r','edgecolor','k','facealpha',1);
    % end

    randomColor = rand(1, 3);
    for i = 1:6
        B = [Fbase(mm, :, i) ;  Ftip(mm, :, i)];
        plot(B(:,1), B(:,2), 'Color', randomColor, 'LineWidth',1.0);
    end
end

daspect([1, 1, 1]);
% title("\bf Filopodia network")
axis off
% axis tight
box on
pp = gcf;
exportgraphics(pp,'Filopodia_Network.png','Resolution',300);


end