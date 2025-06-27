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
t_check = 120;                %% Frequency of checing convergence
MinConvergeTime = 600;        %% minimum expected time of convergence
MinSaveTime = 0;              %% Time to start saving data
FreqSave = 10;                %% How often data is being saved
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
Nc = 20;                      %% Number of cells in a row
MinStableValue = 2;           %% Minimum of changing cells to establish convergence
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
F_mu = 3.50;                 %% (Previous idea): mean length of filopodia
a_apical = 0.1;              %% scaling factor for the apical distribution of Delta
a_basal = 0.1;               %% scaling factor for the apical distribution of Delta

FixVals = [N_mu, D_mu, F_mu, a_apical, a_basal];

%% Parameters for generating filopodia length (prefixed by 'f' to distinguish them from other parameters)
fdelta = 2.7e-3;        %% Half-size of actin monomer
fa0 = 10;               %% G-acting concentration at the leading edge
fKbT = 4.1e-3;          %% Thermal energy
fN = 15;                %% Number of filament in a filopodia bundle
fD = 5;                 %% Effective G-acting diffusion
feta = 20;              %% Geometric conversion coefficient
fFt = 20;               %% Membrane tension
ftau = 2000;            %% 2046.8;     Drag coefficient of filopodia
fK0 = 10;               %% base value for G-acting assembly rate
fm = 0.1;               %% Exponent of retraction rate for K_on
fw = 0.01;              %% exponent of the protrusion rate
fk_sigma = 6e-3;        %% Variation in the membrane force Fm
fNum_filop = Nc*Nc*6;   %% Total number of filopodia in the tissue
fdt = 0.1;              %% Time step for simulating the filopodia
ftf = LifeTime*60;      %% Filopodia life time in seconds
fselectIndx = 60/fdt;   %% Index for selecting filopodia data to be used in the model simulation
fFs = fKbT*fN/fdelta;   %% Membrane tension for which actin polymerization stalls
fP =60;                 %% Period of myosin contraction

LentPar = [fdelta, fa0, fKbT,fN, fD, feta, fFt, ftau, fK0, fm, fw,fk_sigma, fNum_filop, fdt, ...
    ftf,fselectIndx, fFs, LifeTime, fP];

% k_mu = [0.0125, 4, 6];        %% Type A: Coefficients for peak myosin force. ie f0 = k_mu*Ft
% k_mu = [0.0125, 3, 6];        %% Type B: Coefficients for peak myosin force. ie f0 = k_mu*Ft
% k_mu = [0.0125, 4, 12];       %% Type C:  Coefficients for peak myosin force. ie f0 = k_mu*Ft
k_mu = [0.0125, 5, 9];          %% Type D: Coefficients for peak myosin force. ie f0 = k_mu*Ft

alpha = [0.3, 0.65, 0.8; 0.25, 0.6, 0.75; 0.2, 0.75, 0.90];   %% Times when myosin switch phase

Tswitch = [0.25, 0.8]*ftf;      %% Times when filopodia switch phase


%% --- Parameters ranges for sample analysis ------------
SampSize = 10;                        %% Sample size
Np = 10;                              %% Number of partition for the varying parameter

apic = linspace(0.0, 1, Np);          %% Range for apical weight
base = linspace(0.0, 1, Np);          %% Range for basal weight
Dvar = linspace(0.1, 1, Np);          %% Range for initial Delta distribution
Nvar = linspace(0.1, 1, Np);          %% Range for initial Notch distribution

%%%%%==========================================================
%%  ------------ Simulate the various functions below ------------
%%%%%==========================================================
F  = F_sigma*randn(Nc*Nc,6) + F_mu;      %% Distribution of filopodia
Ang = pi*rand(Nc*Nc,6);                  %% Distribution filopodia angles
Notch = N_sigma*randn(Nc*Nc,1) + N_mu;   %% Initial distribution of Notch
Delta = D_sigma*randn(Nc*Nc, 1) + D_mu;  %% Initial distribution of Delta

%%% --- Uncomment the function you would want to run below --------
%%% To run any of the first four functions, you would need to define the filopidia and it angle data 
%%% as well as the Notch and Delta level. You will therefore have to uncomment the above as well .

tic
% [VertCod,CenterCod, Pcell] = TwoDGeom(Nc,r);
% [DirVec, Fbase,FilopDisc, JunctDisc,Ftip] = Filop_vectors(Nc, r, F, Ang, CellsAway);
[DinJ, DinF, FCount] = DeltaIn(Nc, r, F, Ang, CellsAway,a_apical, a_basal, Delta);
% [Lsub] = FilopLent(LentPar, k_mu, alpha, Tswitch);
% [pp, FData]=GenTrendFilopDyn(LentPar,TimeP,ModelP,ContP,N_mu, D_mu,a_apical, a_basal,k_mu, alpha, Tswitch);
Time = toc