clear; close all; clc;
%%%%%%%%%=======================================================================================
% % % global  tf dt t_check MinConvergeTime MinSaveTime FreqSave MaxIt
% % % global R_N R_D mu rho a b h k
% % % global Nc MinStableValue F_rate r N_sigma D_sigma F_sigma N_thresh Rreach
% % % global a_apical a_basal F_mu
%%%%%%%%%=======================================================================================

%%  ------ Time parameters = [tf, dt, t_check, MinConvergeTime, MinSaveTime, FreqSave, MaxIt]-----
tf = 3000;                    %% final simulation time for the model
dt = 1.0;                     %% Temporal resolution (time step)
t_check = 60;                 %% Frequency of checing convergence
MinConvergeTime = 600;        %% minimum expected time of convergence
MinSaveTime = 0;              %% Time to start saving data
FreqSave = 1;                 %% How often data is being saved
MaxIt = tf/dt+1;              %% Number of time iteration
LifeTime = 10;                %% Filopodia lifetime in minutes

TimeP = [tf, dt, t_check,MinConvergeTime, MinSaveTime, FreqSave, MaxIt,LifeTime];

%%  ------Notch-Delta Model parameters = [R_N, R_D, mu, rho, a, b, h, k];--------------------
R_N = 0.1;              %% Rate of production of Notch in the cell
R_D = 0.1;              %% Rate of production of Delta in the cell
mu = 0.002;             %% Rate of decay of Notch in the cell
rho = 0.002;            %% Rate of decay of Delta in the cell
a = 0.01;               %% constant in the differential equation for Notch
b = 100;                %% constant in the differential equation for Delta
h = 9;                  %% Exponent of Notch in the Delta equation
k = 9;                  %% Exponent of Delta in the Notch equation

ModelP = [R_N, R_D, mu, rho, a, b, h, k];

%%  --Control Parameters= [Nc,MinStableValue,F_rate,r, N_sigma, D_sigma, F_sigma, N_thresh,Rreach]; -
Nc = 15;                      %% Number of cells in a row
MinStableValue = 1.6;           %% Minimum of changing cells to establish convergence
F_rate = 0.01;                %% (Previous idea): random rate use to update the filopodia distribution 
r = 1;                        %% Average cell radii
N_sigma =0.01;                %% Standard Deviation of Notch
D_sigma =0.01;                %% Standard Deviation of  Delta
F_sigma = 0.3;                %% (Previous idea): Standard of filopodia 
N_thresh = 1.0;               %% Threshold of Notch to distiguish sop from epithelia cells
CellsAway = 3;                %% Maximum number of cells any filopodia can reach

ContP = [Nc,MinStableValue,F_rate,r, N_sigma, D_sigma, F_sigma, N_thresh,CellsAway];

%% ----- Fixed and most influencing parameters = [N_mu, D_mu, a_apical, a_basal,F_mu] ----------
N_mu = 1.0;                  %% mean density of Notch
D_mu = 1.0;                  %% mean density of  Delta
F_mu = 2.50;                 %% (Previous idea): mean length of filopodia
a_apical = 0.1;              %% scaling factor for the apical distribution of Delta
a_basal = 0.1;               %% scaling factor for the apical distribution of Delta

FixVals = [N_mu, D_mu, F_mu, a_apical, a_basal];

%% ------- Parameter values for simulating filopodia dynamic model ----------
Ft = 20;                  %% Membrane tension
tau = 2400;               %% Filopodia viscosity
delta = 2.7e-3;           %% Half monomer size of actin filament
N = 15;                   %% Number of actin filaments in filopodia
D = 5;                    %% Diffusion coefficient of G-actin
K0 = 10;                  %% Based G-actin assembly rate
KbT = 4.1e-3;             %% Thermal Energy
eta = 20;                 %% Viscosity coefficient
m = 100;                    %% Retraction rate parameter
ftf = LifeTime*60;        %% Filopodia lifetime in seconds
fdt = 0.01;               %% time step for simulating the filopodia model
k_sigma = 3.6;            %% Variability in filopodia length
Num_filop = Nc*Nc*6;      %% Total number of filopodia in the tissue
selectIndx = 60/fdt;      %% Index for selecting filopodia data to be use for model simulation
PauseK0_percent = 0.14;   %% Percentage of G-actin assembly rate during pausing
Fm = 20;                  %% Myosin pulling force
base_G_actin = 20;             %% 20 fixed Basal concentration of G-actin (produced in the cell)

LentPar = [Ft, tau, delta, N,D,K0, KbT, eta, m, ftf, fdt,k_sigma, Num_filop,selectIndx,PauseK0_percent,... 
           Fm, base_G_actin];

%%%%%%%%%============================================================
F  = F_sigma*randn(Nc*Nc,6) + F_mu;      %% Distribution of filopodia
% Ang = pi*rand(Nc*Nc,6);                  %% Distribution filopodia angles
% Ang_NW = ones(Nc*Nc,1)*[(7/6)*pi, pi, (1/2)*pi, (1/6)*pi, 0, (1/3)*pi];
% Ang_NW = ones(Nc*Nc,1)*[(-2/3)*pi, pi, (2/3)*pi, (1/3)*pi, 0, (-1/3)*pi];
% Ang_NE = ones(Nc*Nc,1)*[(2/3)*pi, (1/3)*pi, 0, (5/3)*pi, (4/3)*pi, pi];

Notch = N_sigma*randn(Nc*Nc,1) + N_mu;   %% Initial distribution of Notch
Delta = D_sigma*randn(Nc*Nc, 1) + D_mu;  %% Initial distribution of Delta

%% --- Parameters ranges for sample analysis ------------
SampSize = 10;                        %% Sample size
Np = 10;                              %% Number of partition for the varying parameter

apic = linspace(0.0, 1, Np);          %% Range for apical weight
base = linspace(0.0, 1, Np);          %% Range for basal weight
Dvar = linspace(0.1, 1, Np);          %% Range for initial Delta distribution
Nvar = linspace(0.1, 1, Np);          %% Range for initial Notch distribution
G_actin = linspace(15, 25, Np);

%%%%%==========================================================
%%  ------------ Simulate the various functions below ---------
%%%%%==========================================================
% % [Ang] = FilopAngles(Nc, r, 'C');
% [Lsub] = FilopLent(LentPar,base_G_actin);
[FData]=GenTimeIteration(LentPar,TimeP,ModelP,ContP,FixVals, base_G_actin);
% tic
% [FData]=GenTimeIteration(LentPar,TimeP,ModelP,ContP,FixVals, base_G_actin, Ang);
% [FData]=GenTimeIteration(LentPar,TimeP,ModelP,ContP,FixVals, base_G_actin);
% [Lsub] = FilopLent(LentPar,base_G_actin);
% SwitchMatrix = SwitchingCells(LentPar,TimeP,ModelP,ContP,FixVals, base_G_actin, SampSize);
% Time = toc



% % f = @(gamma) (gamma.*log(gamma) - 2.5271*gamma + 3.8);
% f = @(rho) (100*log(30*rho+28) -60*log(rho) -408.5022);
% lower = fzero(f,[1e-8 , 1.8]);
% upper = fzero(f,[1.8, 20]);
% fprintf('Lower limit = %.4f\n', lower)
% fprintf('Upper limit = %.4f\n', upper)


