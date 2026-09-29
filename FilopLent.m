function [Lsub] = FilopLent(LentPar,G_actin)
% rng(43)
% % LentPar = [Ft, tau, delta, N,D,K0, KbT, eta, m, ftf, fdt,k_sigma, Num_filop,selectIndx,PauseK0_percent];
Ft = LentPar(1);                    %% Membrane tension         
tau = LentPar(2);                   %% Filopodia viscosity
delta = LentPar(3);                 %% Half monomer size of G actin
N = LentPar(4);                     %% Number of actin filaments in filopodia
D = LentPar(5);                     %% Diffusion coefficient of G-actin
K0 = LentPar(6);                    %% Base G-actin assembly rate
KbT = LentPar(7);                  %% Thermal Energy
eta = LentPar(8);                  %% Viscosity coefficient
m = LentPar(9);                    %% Retraction rate parameter
ftf = LentPar(10);                  %% Filopodia lifetime in seconds
fdt = LentPar(11);                  %% time step for simulating the filopodia model
k_sigma = LentPar(12);              %% Variability in filopodia length
Num_filop = LentPar(13);            %% Total number of filopodia in the tissue
selectIndx = LentPar(14);           %% Index for selecting filopodia data to be use for model simulation
PauseK0_percent = LentPar(15);      %% Percentage of G-actin assembly rate during pausing
Fm = LentPar(16);
% G_actin = LentPar(17);

kappa = k_sigma*randn(Num_filop,1); %% variation in the filopodia length

a0 = G_actin;

LifetimeInMinutes = ftf/60;

if LifetimeInMinutes <= 10
    tt0 = 0.50;
    tt1 = 0.50;
else
    tt0 = 5/LifetimeInMinutes;
    tt1 = (LifetimeInMinutes - 5)/LifetimeInMinutes;
end
Tswitch = [tt0, tt1]*ftf;    

t = 0:fdt:ftf;                      %% Simulation time of filopodia
N0 = (Ft*delta)/KbT;                %% Minimum number of filopodia bundle to support filopodia
L0 = zeros(Num_filop, 1);           %% Initial length of filopodia

%%%%%% ----- Defining the protrusion (Vp) and retraction (Vr) rates ----------------
Vp = @(L, Kon, a0)( (Kon*delta*a0)*(1- (Kon*L*N)./(Kon*L*N + D*eta*exp(N0/N) ) )*exp(-(N0/N)) );

Vr = (Fm + Ft)/tau;

%%%%% -----Initialing vecors to store data ---------
Lmat = zeros(Num_filop, length(t));

%%%%%% ------------ Solving the filopodia model ------------------------
for n = 2: length(t)

    Kon = KonBar(t(n), Tswitch, K0, m,PauseK0_percent);     %% G-actin assemble rate

    Lnew = L0 + fdt*( (Vp(L0, Kon, a0) - Vr) );

    Lmat(:,n) = Lnew + kappa;
    
    L0 = Lnew;
    
end

Scale_Lmat = Lmat/5;                     % scaling the length to into cell radii
Lsub = Scale_Lmat(:,1:selectIndx:end);   % Extacting filopodia data for the main simulation 


%%%%%%%=====================================================================
%%%%%%%% --------------- Visualizing the filopodia dynamics ----------------
% figure
% figure('Position',[100 100 600 380])
% 
% Scale_err = std(Scale_Lmat);
% Scale_md = mean(Scale_Lmat);
% maxLent = max(Scale_md);
% 
% errorbar( t, Scale_md,  Scale_err, 'c')
% hold on
% plot(t, Scale_md, 'm', 'LineWidth',2);
% yline(maxLent, 'b --','LineWidth',2)
% 
% txt = ['G-actin ', '$(a_0) = ', num2str(a0), '$', newline...,
%        'Lifetime = ', num2str(ftf/60),' min', newline...,
%        'Avg. length = ', num2str(maxLent)];
% 
% text(0.3,0.150,txt,'Units','normalized','Interpreter','latex','FontSize',16,'FontWeight','bold', ...
%     'BackgroundColor','w','EdgeColor','k')
% 
% ylim([0,3.0])
% xlim([0,ftf])
% xlabel('\bf Time [sec]', 'FontSize',14, 'FontWeight','bold'); 
% ylabel('\bf Length [Cell radii]', 'FontSize',14, 'FontWeight','bold')
% pb = gcf;
% exportgraphics(pb,'Filop_Length.png','Resolution',300);
% 



%%%%%%%%%=========================================================================
    function Kon = KonBar(t, Tswitch, K0, m, PauseK0_percent)

        t0 = Tswitch(1);                 %% Time when filopodia start pausing
        t1 = Tswitch(2);                 %% Time when filopodia start retracting

        if t <= t0                       %% G-actin assembly during protrusion
            Kon = K0;

        elseif (t >t0 && t< t1)          %% G-actin assembly during pausing
            Kon = PauseK0_percent*K0;

        else                             %% G-actin assembly during retraction
            Kon = K0*exp(-m*(t - t1));
        end

    end


end
