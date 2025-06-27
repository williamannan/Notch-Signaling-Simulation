function [pp, FData]=GenTrendFilopDyn(LentPar,TimeP,ModelP,ContP,N_mu, D_mu,a_apical, a_basal,k_mu, alpha, Tswitch)

%% ---------Defining model parameters ----------------------
%%% TimeP = [tf, dt, t_check,MinConvergeTime, MinSaveTime, FreqSave, MaxIt];
% tf = TimeP(1);               %% Final Simulation time      
dt = TimeP(2);                 %% Time step
t_check = TimeP(3);            %% Time to check convergence      
MinConvergeTime = TimeP(4);    %% Minimum time to start saving simulated data
MinSaveTime = TimeP(5);        %% Minimum time the pattern may stabilize
FreqSave = TimeP(6);           %% How often simulated data is saved
MaxIt = TimeP(7);              %% NUmber of time iteractions
LifeTime = TimeP(8);           %% Filopodia life time in minutes 

%%% ModelP = [R_N, R_D, mu, rho, a, b, h, k];
R_N = ModelP(1);              %% Rate of Notch prodaction
R_D = ModelP(2);              %% Rate of Delta prodaction
mu = ModelP(3);               %% Rate of decay of Notch
rho = ModelP(4);              %% Rate of decay of Delta
a = ModelP(5);                %% constant in the Notch equation
b = ModelP(6);                %% constant in the Delta equation
h = ModelP(7);                %% constant in the Delta equation
k = ModelP(8);                %% constant in the Notch equation

%%% ContP = [Nc,MinStableValue,F_rate,r, N_sigma, D_sigma, F_sigma, N_thresh,CellsAway];
Nc = ContP(1);                %% Number of cells in a row
MinStableValue = ContP(2);    %% Minimum number of changing cells to establish convergence
% F_rate = ContP(3);          %% Not useful when using filopodia dynamics model
r = ContP(4);                 %% Cell radius 
N_sigma = ContP(5);           %% Standard deviation for Notch
D_sigma = ContP(6);           %% Standard deviation for Delta
% F_sigma = ContP(7);         %% Not useful when using filopodia dynamics model
N_thresh = ContP(8);          %% Notch threshold to distinguish SOP from epithelia 
CellsAway = ContP(9);         %% Maximum # of cell diameters a filopodium can reach

%%% ------- Geometry of the cells ---------
FLent = FilopLent(LentPar,k_mu, alpha, Tswitch);     %% length distribution of filopodia
LifeCount = 1;                                     %% keeping track of the number of filopodia life times      

% F  = F_sigma*randn(Nc*Nc,6) + F_mu;              %% Not useful when using filopodia dynamics model

Ang = pi*rand(Nc*Nc, 6);                           %% Distribution filopodia angles
Notch = N_sigma*randn(1, Nc*Nc) + N_mu;            %% Initial distribution of Notch
Delta = D_sigma*randn(1, Nc*Nc,1) + D_mu;          %% Initial distribution of Delta

%%%%%%%=====================================================================
%% ---------  Main time stepping loop to update Delta and Notch  ------------
%%%%%%%=====================================================================
CellState = zeros(MaxIt, Nc*Nc);         %% Initializing state vector for the cell
ScNotch = zeros(MaxIt, Nc*Nc + 1);       %% Initialized matrix for saving scaled Notch
ANotch = zeros(MaxIt, Nc*Nc);            %% Initialized matrix for saving actual Notch

cntSave = 1;                             %% counter to keep track of saved data
Conv_count = 1;                          %% counter to keep track of convergence test
t0 = 0;                                  %% times to stored data

MeanDin = zeros(MaxIt, 3);               %% To keep track of mean Din at each time step                     
UnstCellCount = zeros(MaxIt,2);          %% To keep track of number of cell changing states

for step = 1: MaxIt
    LocInd = step + 1 - (LifeCount - 1)*LifeTime;   %% Index for selecting filopodia data to be used

    F = reshape(FLent(:,LocInd), Nc*Nc, 6);  %% Organizing the filopodia data 
    
    %%%%%%%************************************************************************
    %%% ---- Computing effective Delta express by neighboring cells --------
    [DinJ, DinF, ~] = ParDeltaIn(Nc, r, F, Ang, CellsAway,a_apical, a_basal, Delta);

    Din = DinJ + DinF;

    MeanDin(step, :) =[step, mean(DinJ), mean(DinF)]; 
    
    %%%%%%%************************************************************************
    %%  --- Updating Notch and Delta distribution at each time step ----
    Notch = Notch + dt*( R_N* (Din.^k)./(a + Din.^k) - mu*Notch);
    Delta = Delta + dt*( R_D* 1./(1 + b*Notch.^h) - rho*Delta);

    %% ---Previous idea on how to update Filopodia distribution at each time step ----
    % F_updateIndex = 100*rand(Nc*Nc, 6) < F_rate*1000;         %% indices for filopodia that need update
    % Fnew = F_sigma*randn(Nc*Nc, 6) + F_mu;                    %% Generating new filopodia for all cells
    % Ang = pi*rand(Nc*Nc, 6);
    % F(F_updateIndex ==1) = Fnew(F_updateIndex ==1);         %% Updating filopodia
    % Ang(F_updateIndex ==1) = NewAng(F_updateIndex ==1);      %% Updating angles
    
    %%%%%%%************************************************************************
    %% ---Current idea on how to update Filopodia distribution at each time step ----
    %%% NOTE: Whenever the life time of the filopodia lapses, we sample new angles and solve the filopodia model
    if rem(step, LifeTime) == 0            
        Ang = pi*rand(Nc*Nc, 6);                       %% sample new filopodia angles
        FLent = FilopLent(LentPar,k_mu, alpha, Tswitch); %% Solve the filopodia model 
        LifeCount = LifeCount + 1;                     %% increase counter for the life time
    end

   %%%%%%%************************************************************************
    %% --- scaling Notch by moving averages and saving data for analysis -------
    if ((step >= MinSaveTime) && (rem(step, FreqSave) == 0) )

        ANotch(cntSave, :) = Notch;      

        if cntSave == 1
            ScNotch(cntSave, :) = [t0, Notch./Notch];
        else
            ScNotch(cntSave, :) = [t0, Notch./mean(ANotch(1:cntSave, :))];
        end

        cntSave = cntSave + 1;
        t0 = t0 + FreqSave;
    end

    %%%%%%%*********************************************************************
    %% ------ Checking for convergence -------------------------
    if ((step >= MinConvergeTime) && rem(step,t_check) == 0)

        StateNotch =  ScNotch(cntSave-1, 2:end);

        CellState(Conv_count, :) = StateNotch <= N_thresh;  %% Binary states of the cells

        %%% computing the number of cells that are changing states over
        if Conv_count > 1
            numOfUnstCells = sum( abs( CellState(Conv_count,:) - CellState(Conv_count-1,:) ) );
            UnstCellCount(Conv_count-1, :) = [step, numOfUnstCells];
        end

        %%%  breaking the time iteration loop if the system is stable
        if ( (Conv_count > 1) && ( numOfUnstCells/(Nc*Nc) <  MinStableValue/100) )
            break
        else
            Conv_count = Conv_count + 1;
        end
        % % % count = count + 1;
    end

end

%%%%%%%************************************************************************
%% ------Output data after one complete time iteration ---------------------
FinalTime = step;                    %% Time taken for stable pattern formation in minutes 
ScNotch = ScNotch(1:cntSave-1, :);   %% Scaled Notch 
ANotch = ANotch(1:cntSave-1, :);     %% Actual Notch 

MeanDin = MeanDin(1:FinalTime, :);   %% mean Din Data
FinalNotch = ScNotch(end,2:end);     %% Scaled Notch distribution at the final time


%%%%%%%%%%=========================================================================
%%%%%%%%%% ============ ---- Visualing output  -------=============================

%% ----- Computing Relative and Local sensory area ------------------------
[~, ~,FilopDisc, ~,~] = Filop_vectors(Nc, r, F, Ang, CellsAway);

FinalSOPIndx = find(FinalNotch <= N_thresh);    %% extracting the indices of SOP cells

LSA = zeros(length(FinalSOPIndx),1);             %% Local sensory area

TissueArea = Nc*Nc*(pi*5^2);                     %% Tissue are radius of cell is roughly 5 micrometers 
RSA = (length(FinalSOPIndx)*28.45)/(TissueArea); %% Relative density

for jj = 1:length(FinalSOPIndx)
    SOPs_in_disc = intersect(FilopDisc{FinalSOPIndx(jj)}, FinalSOPIndx);
    DiscArea = length(FilopDisc{FinalSOPIndx(jj)})*(pi*5^2);

    LSA(jj) = (length(SOPs_in_disc)*28.45)/DiscArea;
end

FData = [FinalTime , RSA, mean(LSA)]; 

%%%%%%%************************************************************************
%%% --- making a boxplot for Local sensory area --------
figure
boxplot(LSA);
str={'SOP'};    
ylabel("\bf LSA")
title(['\bf [Time, RSA, LSA] = ',num2str(FData)])
set(gca,'XTickLabel',str)
pp(1) = gcf;
exportgraphics(pp(1),'boxplot_LSA.png','Resolution',300);

%%%%%%%******************************************************************
%%% ---------- Plotting the trend of mean Din --------------------------------
figure
plot(MeanDin(:,1),MeanDin(:,2), 'r',MeanDin(:,1), MeanDin(:,3), 'b', 'LineWidth',1.5 )
title('\bf D_{in} profile ')
xlabel("\bf Time")
ylabel("\bf Din")
legend("Junctional", "Basal",'Location','east')
axis tight
box on
pp(3) = gcf;
exportgraphics(pp(3),'Din_profile.png','Resolution',300);


%%%%%%%******************************************************************
%%% ---------- Plotting the Notch profile--------------------------------
sopInd = FinalNotch <= N_thresh;     %% index for the SOP cells
epInd = FinalNotch > N_thresh;       %% Index for epithelia cells

sopNotch =  ANotch(:, sopInd);       %% Extracting actual Notch for SOP cells
epNotch = ANotch(:, epInd);          %% Extracting actual Notch for epithelia cells

sopNotch = sopNotch./max(sopNotch);  %% Scaling the SOP Notch the maximum
epNotch = epNotch./max(epNotch);     %% Scaling the epithelia Notch the maximum

MeanSOPNotch = mean(sopNotch, 2);    %% Taking average SOP Notch accross the tissue
MeanEPNotch = mean(epNotch, 2);      %% Taking average epithelia Notch accross the tissue

Tt = ScNotch(:,1);                   %% Simulation time

figure
plot(Tt, MeanSOPNotch,'r', Tt, MeanEPNotch, 'b', 'LineWidth',1.5);
title('\bf Notch profile')
xlabel("\bf Time")
ylabel("\bf Notch")
legend("SOP", "Epi",'Location','east')
axis tight
box on
pp(4) = gcf;
exportgraphics(pp(4),'Notch_profile.png','Resolution',300);

%%%%%%%******************************************************************
%%% ---------- Plotting the final pattern formed ------------------------
[~,~, Pcell] = TwoDGeom(Nc,r);
figure
plot(Pcell,'FaceColor','c','edgecolor','k','facealpha',1);
hold on
for m = 1:Nc*Nc
    if FinalNotch(m) < N_thresh
        plot(Pcell(m),'FaceColor','r','edgecolor','k','facealpha',1);
    end
end
title('\bf Final pattern')
daspect([1, 1, 1])
axis tight
box on
pp(5) = gcf;
exportgraphics(pp(5),'2DHexPattern.png','Resolution',300);


end