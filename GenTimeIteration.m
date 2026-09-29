function [FData]=GenTimeIteration(LentPar,TimeP,ModelP,ContP,FixVals, G_actin)
% rng(43)
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

%%%=------FixVals = [N_mu, D_mu, F_mu, a_apical, a_basal];
N_mu = FixVals(1);
D_mu = FixVals(2);
% F_mu = FixVals(3);
a_apical = FixVals(4);
a_basal =FixVals(5);

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

%%% ------- Distribution of filopodia length ---------
FLent = FilopLent(LentPar, G_actin);            %% length distribution of filopodia
LifeCount = 1;                                        %% keeping track of the number of filopodia life times

Ang = pi*rand(Nc*Nc, 6);                              %% Distribution filopodia angles
Notch = N_sigma*randn(1, Nc*Nc) + N_mu;               %% Initial distribution of Notch
Delta = D_sigma*randn(1, Nc*Nc,1) + D_mu;             %% Initial distribution of Delta

%%%%%%%=====================================================================
%% ---------  Main time stepping loop to update Delta and Notch  -----------
%%%%%%%=====================================================================
ScNotch = zeros(MaxIt, Nc*Nc + 1);       %% Initialized matrix for saving scaled Notch
ANotch = zeros(MaxIt, Nc*Nc);            %% Initialized matrix for saving actual Notch

ApicDinMat = zeros(MaxIt, Nc*Nc);
FilopDinMat = zeros(MaxIt, Nc*Nc);
AllTimes = zeros(MaxIt, 1);

cntSave = 1;                             %% counter to keep track of saved data
Conv_count = 1;                          %% counter to keep track of convergence test
t0 = 0;                                  %% times to stored data

MeanDin = zeros(MaxIt, 3);               %% To keep track of mean Din at each time step

UnstCellCount = zeros(MaxIt,2);          %% To keep track of number of cell changing states

PrevCellState = [];
EverChanged = false(1,Nc*Nc);

for step = 1: MaxIt
    LocInd = step + 1 - (LifeCount - 1)*LifeTime;   %% Index for selecting filopodia data to be used

    F = reshape(FLent(:,LocInd), Nc*Nc, 6);  %% Organizing the filopodia data

    %%%%%%%************************************************************************
    %%% ---- Computing effective Delta express by neighboring cells --------
    [DinJ, DinF, ~] = ParDeltaIn(Nc, r, F, Ang, CellsAway,a_apical, a_basal, Delta);

    ApicDinMat(step, :) = DinJ;
    FilopDinMat(step, :) = DinF;
    AllTimes(step) = step;

    Din = DinJ+ DinF;

    MeanDin(step, :) =[step, mean(DinJ), mean(DinF)];

    %%%%%%%************************************************************************
    %%  --- Updating Notch and Delta distribution at each time step ----
    Notch = Notch + dt*( R_N* (Din.^k)./(a + Din.^k) - mu*Notch);
    Delta = Delta + dt*( R_D* 1./(1 + b*Notch.^h) - rho*Delta);

    %%%%%%%********************************************************
    %%% ---updating Filopodia distribution at each time step ----
    %%% NOTE: Whenever the life time of the filopodia lapses, 
    %%% we sample new angles and solve the filopodia model

    if rem(step, LifeTime) == 0
        Ang = pi*rand(Nc*Nc, 6);                       %% sample new filopodia angles
        FLent = FilopLent(LentPar, G_actin);
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

    %%%%%%%*************************************
    %% ------ Checking for convergence ---------
    %%% This section of the code has been modified.
    if (step >= MinConvergeTime)

        StateNotch = ScNotch(cntSave-1,2:end);

        CurrentCellState = StateNotch <= N_thresh;

        %%% First time after MinConvergeTime: just initialize
        if isempty(PrevCellState)

            PrevCellState = CurrentCellState;

        else

            %%% Individual cells that changed since previous time step
            Changed = abs(CurrentCellState - PrevCellState);

            %%% Keep track of every cell that has changed during this interval
            EverChanged = EverChanged | Changed;  %% Using logical OR statement 

            %%% Update previous state
            PrevCellState = CurrentCellState;

        end

        %%% Every t_check iterations, evaluate convergence
        if rem(step-MinConvergeTime,t_check) == 0 && step > MinConvergeTime

            numOfUnstCells = sum(EverChanged);

            UnstCellCount(Conv_count,:) = [step,numOfUnstCells];

            disp([step,numOfUnstCells, Nc^2]) 

            if numOfUnstCells/(Nc*Nc) < MinStableValue/100
                break
            end

            Conv_count = Conv_count + 1;

            % Start a fresh interval
            EverChanged = false(1,Nc*Nc);

        end

    end

    %%%%%%%*************************************
    % %% ---- Checking for convergence ---------
    % %%% This was the previous way of checking convergence.

    % if ((step >= MinConvergeTime) && rem(step,t_check) == 0)
    %
    %     StateNotch =  ScNotch(cntSave-1, 2:end);
    %
    %     CurrentCellState = StateNotch <= N_thresh;
    %     CellState(Conv_count, :) = CurrentCellState ;%StateNotch <= N_thresh;  %% Binary states of the cells
    %
    %     disp([step, sum(CurrentCellState), Nc^2])
    %
    %     %%% computing the number of cells that are changing states over
    %     if Conv_count > 1
    %         numOfUnstCells = sum( abs( CellState(Conv_count,:) - CellState(Conv_count-1,:) ) );
    %         UnstCellCount(Conv_count-1, :) = [step, numOfUnstCells];
    %     end
    %
    %     %%%  breaking the time iteration loop if the system is stable
    %     if ( (Conv_count > 1) && ( numOfUnstCells/(Nc*Nc) <  MinStableValue/100) )
    %         break
    %     else
    %         Conv_count = Conv_count + 1;
    %     end
    %     % % % count = count + 1;
    % end

end

%%%%%%%************************************************************************
%% ------Output data after one complete time iteration ---------------------
FinalTime = step;                     %% Time taken for stable pattern formation in minutes
ScNotch = ScNotch(1:cntSave-1, :);   %% Scaled Notch
ANotch = ANotch(1:cntSave-1, :);     %% Actual Notch

MeanDin = MeanDin(1:FinalTime, :);   %% mean Din Data
FinalNotch = ScNotch(end,2:end);     %% Scaled Notch distribution at the final time

ApicDinMat = ApicDinMat(1:FinalTime, :);  %% Din from apical interaction
FilopDinMat= FilopDinMat(1:FinalTime, :); %% Din from filopodia interaction
AllTimes = AllTimes(1: FinalTime);        %% Corresponding times for Din    

UnstCellCount = UnstCellCount(1: Conv_count, :);  %% Number of switching cells


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
% figure
figure('Position',[100 100 600 400])
boxplot(LSA);
txt = ['G-actin ', '$(a_0) = ', num2str(G_actin), '$', newline...,
    'Lifetime = ', num2str(LentPar(10)/60),' min' newline...,
    'Avg. LSA = ', num2str(mean(LSA)), newline...,
    'RSA = ', num2str(RSA)];

text(0.02,0.80,txt,'Units','normalized','Interpreter','latex','FontSize',12, ...
    'FontWeight','bold', 'BackgroundColor','w','EdgeColor','k')

str={'SOP'};
ylabel("\bf LSA")
% title(['\bf [Time, RSA, LSA] = ',num2str(FData)])
set(gca,'XTickLabel',str)
pp(1) = gcf;
exportgraphics(pp(1),'boxplot_LSA.png','Resolution',300);

%%%%%%%******************************************************************
SOPInd = FinalNotch <= N_thresh;     %% index for the SOP cells
EpInd = FinalNotch > N_thresh;       %% Index for epithelia cells

Num_sop_cells = sum(SOPInd);
Num_epit_cells = sum(EpInd);

%%% --- Plotting the average Din expressed by each junctional neighbor ---
MeanApicSOPDin = mean(ApicDinMat(:,SOPInd),2)/Num_sop_cells;
MeanApicEpithDin = mean(ApicDinMat(:,EpInd),2)/Num_epit_cells;

mark_spacing = ceil(length(AllTimes)/20);

% figure
figure('Position',[100 100 600 400])
plot(AllTimes, MeanApicSOPDin, 'm *-', 'LineWidth', 2.0, 'MarkerIndices',...
    1:mark_spacing:length(AllTimes), 'DisplayName', 'Apical SOP D_{in}');
hold on
plot(AllTimes, MeanApicEpithDin, 'b o--', 'LineWidth', 2.0, 'MarkerIndices', ...
    1:mark_spacing:length(AllTimes), 'DisplayName', 'Apical Epithelia D_{in}');


xlabel("\bf Time [min]", 'FontSize',13)
ylabel("\bf Average Junctional D_{in} ", 'FontSize',13)
% title('\bf Junctional: Average D_{in} affecting each cell', 'FontSize',13)
legend('show', 'Location','west')
axis tight
box on
pp(2) = gcf;
exportgraphics(pp(2),'Apical_Din_Per_Neighbor.png','Resolution',300);

%%%%%%%*******************************************
%%% ----- Plotting the average Din expressed by each filopodia neighbor -------
MeanFilopSOPDin = mean(FilopDinMat(:,SOPInd), 2)/Num_sop_cells;
MeanFilopEpithDin = mean(FilopDinMat(:,EpInd),2)/Num_epit_cells;

% figure
figure('Position',[100 100 600 400])

plot(AllTimes, MeanFilopSOPDin, 'm *-', 'LineWidth', 2.0, 'MarkerIndices',...
    1:mark_spacing:length(AllTimes), 'DisplayName', 'Basal SOP D_{in}');
hold on
plot(AllTimes, MeanFilopEpithDin, 'b o--', 'LineWidth', 2.0, 'MarkerIndices', ...
    1:mark_spacing:length(AllTimes), 'DisplayName', 'Basal Epithelia D_{in}');

xlabel("\bf Time [min]", 'FontSize',13)
ylabel("\bf Average Filopodia D_{in} ", 'FontSize',13)
% title('\bf Filopodia: Average D_{in} affecting each cell', 'FontSize',13)
legend('show', 'Location','northwest')
axis tight
box on
pp(3) = gcf;
exportgraphics(pp(3),'Basal_Din_Per_Neighbor.png','Resolution',300);

%%%%%%%**************************************
%%% ----- Plotting the Notch profile---------
sopNotch =  ANotch(:, SOPInd);       %% Extracting actual Notch for SOP cells
epNotch = ANotch(:, EpInd);          %% Extracting actual Notch for epithelia cells

sopNotch = sopNotch./max(sopNotch);  %% Scaling the SOP Notch the maximum
epNotch = epNotch./max(epNotch);     %% Scaling the epithelia Notch the maximum

MeanSOPNotch = mean(sopNotch, 2);    %% Taking average SOP Notch accross the tissue
MeanEPNotch = mean(epNotch, 2);      %% Taking average epithelia Notch accross the tissue

Tt = ScNotch(:,1);                   %% Simulation time

% figure
figure('Position',[100 100 600 400])
plot(Tt, MeanSOPNotch,'b --', Tt, MeanEPNotch, 'm -', 'LineWidth',2.0);

% title('\bf Average Notch level', 'FontSize', 13)
xlabel("\bf Time [min]", 'FontSize', 13)
ylabel("\bf Average Notch ", 'FontSize', 13)

legend("SOP", "Epi",'Location','north')
axis tight
box on
pp(4) = gcf;
exportgraphics(pp(4),'Notch_profile.png','Resolution',300);

%%%%%%%*******************************************
%%% --- Ploting the switching cell profile -------
figure('Position',[100 100 600 400])

plot(UnstCellCount(:,1), UnstCellCount(:,2),'b o--', 'LineWidth',2.0);

% title('\bf Average Notch level', 'FontSize', 13)
xlabel("\bf Time [min]", 'FontSize', 13)
ylabel("\bf Num of switching cells ", 'FontSize', 13)

axis tight
box on
pp(5) = gcf;
exportgraphics(pp(5),'Switching_Cells.png','Resolution',300);


%%%%%%%*******************************************
%%% ---- Plotting the final pattern formed -------
[~,~, Pcell] = TwoDGeom(Nc,r);
% figure
figure('Position',[100 100 600 400])
plot(Pcell,'FaceColor','c','edgecolor','k','facealpha',1);
hold on
for m = 1:Nc*Nc
    if FinalNotch(m) < N_thresh
        plot(Pcell(m),'FaceColor','r','edgecolor','k','facealpha',1);
    end
end
% title('\bf Final pattern')
daspect([1, 1, 1])
axis off
box on
pp(6) = gcf;
exportgraphics(pp(6),'2DHexPattern.png','Resolution',300);

end


