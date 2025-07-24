function [Lsub] = FilopLent(LentPar, alpha,Tswitch)

%%% LentPar = [Ft, kf, tau, delta, a0, r_filop, N,D,K0, kbT, eta, m,alpha0, alpha1, ftf, fdt];
Ft = LentPar(1);                    %% Membrane tension         
kf = LentPar(2);                    %% Filopodia stiffness
tau = LentPar(3);                   %% Filopodia viscosity
delta = LentPar(4);                 %% Half monomer size of G actin
a0 = LentPar(5);                    %% fixed Basal concentration of G-actin (produced in the cell)
r_filop = LentPar(6);               %% Scaling factor for myosin force Fm (assumed)
N = LentPar(7);                     %% Number of actin filaments in filopodia
D = LentPar(8);                     %% Diffusion coefficient of G-actin
K0 = LentPar(9);                    %% Base G-actin assembly rate
KbT = LentPar(10);                  %% Thermal Energy
eta = LentPar(11);                  %% Viscosity coefficient
m = LentPar(12);                    %% Retraction rate parameter
ftf = LentPar(13);                  %% Filopodia lifetime in seconds
fdt = LentPar(14);                  %% time step for simulating the filopodia model
k_sigma = LentPar(15);              %% Variability in filopodia length
Num_filop = LentPar(16);            %% Total number of filopodia in the tissue
selectIndx = LentPar(17);           %% Index for selecting filopodia data to be use for model simulation
PauseK0_percent = LentPar(18);      %% Percentage of G-actin assembly rate during pausing

kappa = k_sigma*randn(Num_filop,1); %% variation in the filopodia length


alpha0 = alpha(1);                  %% Base myosin pulling force
alpha1 = alpha(2);                  %% Coefficient of local myosin force due to filopodia length

t = 0:fdt:ftf;                      %% Simulation time of filopodia
N0 = (Ft*delta)/KbT;                %% Minimum number of filopodia bundle to support filopodia
L0 = 0;                             %% Initial length of filopodia
x0 = 0;                             %% Initial retrograde displacement

%%%%%% ----- Defining the protrusion rate Vp and the myosin pulling force Fm   ----------------
Vp = @(L, Kon, a0bar)( (Kon*delta*a0bar)*(1- (Kon*L*N)./(Kon*L*N + D*eta*exp(N0/N) ) )*exp(-(N0/N)) );
Fm = @(L)(alpha0 + alpha1*L);

%%%%% ---------------Initialing vecors to store data ---------------------------
Lmat = zeros(Num_filop, length(t));
xMat = zeros(1, length(t));
VrMat = zeros(1, length(t));
KonMat = zeros(1, length(t));
VpMat = zeros(1, length(t));

%%%%%% ------------ Solving the filopodia model ------------------------
for n = 2: length(t)
    Kon = KonBar(t(n), Tswitch, K0, m,PauseK0_percent);     %% G-actin assemble rate
    a0bar = a0 + r_filop*Fm(L0);                             % a0 + r_filop*Fm(L0)   
    Vr = (Fm(L0) + Ft - kf*x0)/tau;
    xnew = x0 + fdt*Vr; 

    Lnew = L0 + fdt*( Vp(L0, Kon, a0bar) - Vr);

    
    Lmat(:,n) = Lnew + kappa;

    xMat(n) = xnew;
    VrMat(n) = Vr;
    KonMat(n) = Kon;
    VpMat(n) = Vp(L0, Kon, a0bar);

    L0 = Lnew;
    x0 = xnew;
end

Scale_Lmat = Lmat/5;                     % scaling the length to into cell radii
Lsub = Scale_Lmat(:,1:selectIndx:end);   % Extacting filopodia data for the main simulation 

% % % %%%%% --------------- Visualizing the filopodia dynamics ----------------
% % % figure
% % % Scale_err = std(Scale_Lmat);
% % % Scale_md = mean(Scale_Lmat);
% % % errorbar( t, Scale_md,  Scale_err, 'm')
% % % hold on
% % % plot(t, Scale_md, 'b', 'LineWidth',2);
% % % ylim([0,2.5])
% % % title('\bf Filopodia dynamics; Lifetime')
% % % xlabel('\bf Time [sec]'); 
% % % ylabel('\bf Filopodia length in cell radii')
% % % pb = gcf;
% % % exportgraphics(pb,'Scale_Filop_Dyn.png','Resolution',300);
% % % 
% % % figure
% % % err = std(Lmat);
% % % md = mean(Lmat);
% % % errorbar( t, md,  err, 'c')
% % % hold on
% % % plot(t, md, 'r', 'LineWidth',2);
% % % ylim([0,12])
% % % title('\bf Actual filopodia length dynamics')
% % % xlabel('\bf Time [sec]'); 
% % % ylabel('\bf Filopodia length in \mu m')
% % % pb = gcf;
% % % exportgraphics(pb,'Actual Filop_Dyn.png','Resolution',300);
% % % 
% % % figure
% % % plot(t(2:end),KonMat(2:end), 'b--','LineWidth',2)
% % % xlabel('\bf Time', 'FontSize',14);
% % % ylabel('\bf G-actin assembly rate K_{on}', 'FontSize',14)
% % % 
% % % figure
% % % plot(t(2:end),-xMat(2:end), 'b--','LineWidth',2)
% % % xlabel('\bf Time', 'FontSize',14);
% % % ylabel('\bf Retrograde displacement x(t)', 'FontSize',14)
% % % 
% % % figure
% % % plot(t(2:end), VpMat(2:end), 'b:',t(2:end), VrMat(2:end), 'm--','LineWidth',2)
% % % xlabel('\bf Time', 'FontSize',14);
% % % ylabel('\bf Protrusion and retraction rates', 'FontSize',14)
% % % legend('Protrusion rate V_p', "Retraction rate V_r")

%%%%%%%%%=========================================================================
    function Kon = KonBar(t, Tswitch, K0, m, PauseK0_percent)

        t0 = Tswitch(1);                 %% Time when filopodia start pausing
        t1 = Tswitch(2);                 %% Time when filopodia start retracting

        if t <= t0                       %% G-actin assembly during protrusion
            Kon = K0;

        elseif (t >t0 && t< t1)          %% G-actin assembly during pausing
            Kon = PauseK0_percent*K0;

        else                             %% G-actin assembly during retraction
            Kon = K0*exp(-m*(t - t0));
        end

    end


end
