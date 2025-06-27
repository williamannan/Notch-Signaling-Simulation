function [Lsub] = FilopLent(LentPar, k_mu, alpha, Tswitch)

delta = LentPar(1);      %% Half-size of actin monomer
a0 = LentPar(2);         %% G-acting concentration at the leading edge
KbT = LentPar(3);        %% Thermal energy
N = LentPar(4);          %% Number of filament in a filopodia bundle
D = LentPar(5);          %% Effective G-acting diffusion
eta = LentPar(6);        %% Geometric conversion coefficient
Ft = LentPar(7);         %% Membrane tension
tau = LentPar(8);        %% Drag coefficient of filopodia
K0 = LentPar(9);         %% base value for G-acting assembly rate
m = LentPar(10);         %% Exponent of retraction rate for K_on
w = LentPar(11);         %% exponent of the protrusion rate
k_sigma = LentPar(12);   %% Variation in the membrane force Fm
Num_filop = LentPar(13); %% Total number of filopodia in the tissue
dt = LentPar(14);        %% Time step for simulating the filopodia
tf = LentPar(15);        %% Filopodia life time in seconds
selectIndx = LentPar(16); %% Selected index for filopodia data
Fs = LentPar(17);        %% Membrane tension for which actin polymerization stalls
LifeTime = LentPar(18);  %% Filopodia life time in minutes
P = LentPar(19);         %% Period of myosin contraction cycle

t = 0:dt:tf;             %% Simulation time

N0 = Ft*delta/KbT;       %% Number of filament to support filopodium
c0 = (1 - (Ft/Fs))^w;
kappa = k_sigma*randn(Num_filop,1);  %% variation 

Vp = @(L, Kon)( (Kon*delta*a0*c0)*(1- (Kon*L*N)./(Kon*L*N + D*eta*exp(N0/N)) ) );


Lmat = zeros(Num_filop, length(t));
Lmat(:,1) = 0;
FmData = zeros(length(t), 1);
MyosinData = zeros(length(t), 1);

for n = 1: length(t)-1
    if t(n) <= Tswitch(1)   % Myosin force and retrograte flow rate during elongation phase
        f0 = Ft*k_mu(1);
        Fm = OneMyosin(t(n), alpha(1, :), f0, P);
        Vr = (Ft + Fm)/tau + kappa;
        FmData(n) = Fm;

    elseif (t(n) > Tswitch(1)) && (t(n) <= Tswitch(2))   % Myosin force and retrograte flow rate during stationary phase
        f0 = Ft*k_mu(2);
        Fm = OneMyosin(t(n), alpha(2, :), f0, P);
        Vr = (Ft + Fm)/tau + kappa;
        FmData(n) = Fm;

    else
        f0 = Ft*k_mu(3);       % Myosin force and retrograte flow rate during retraction phase
        Fm = OneMyosin(t(n), alpha(3, :), f0, P);
        Vr = (Ft + Fm)/tau + kappa;
        FmData(n) = Fm;
    end

    
    [Kon] = KonBar(t(n),Tswitch,K0, m);

    Lnew =  Lmat(:,n) + dt*( Vp(Lmat(:,n), Kon) - Vr );
    Lnew(Lnew<0) = 0;            % setting the negative length to zero
    Lmat(:,n+1) =  Lnew;

    MyosinData(n) = OneMyosin(t(n), alpha(1, :), k_mu(1)*Ft, P);
end

Lmat = Lmat/5;                    % scaling the length to into cell radii
Lsub = Lmat(:,1:selectIndx:end);  % Extacting filopodia data for the main simulation 

% figure
% plot(Lmat)

figure
err = std(Lmat);
md = mean(Lmat);
errorbar( t, md,  err, 'c')
hold on
plot(t, md, 'r', 'LineWidth',2);
title(['\bf Filopodia dynamics; Lifetime = ', num2str(LifeTime),' minutes'])
xlabel('\bf Time [sec]'); 
ylabel('\bf Filopodia length in cell radii')
pb = gcf;
exportgraphics(pb,'Filop_Dyn.png','Resolution',300);

figure
plot(t, FmData, 'r', 'LineWidth',1.5);
title('\bf Myosin force ')
xlabel('\bf Time [sec]'); ylabel('\bf Fm')
pf = gcf;
exportgraphics(pf,'Myosin_Force.png','Resolution',300);

figure
plot(t, MyosinData, 'r', 'LineWidth',1.5);
title('\bf Myosin Profile ')
ylim([0, 0.5])
xlabel('\bf Time [sec]');
ylabel('\bf Myosin force F_m')
box on; 
pm = gcf;
exportgraphics(pm,'Myosin.png','Resolution',300);



%%%%%%------- User define functions for the simulation ---------
function [Kon] = KonBar(t,Tswitch,K0, m)
        t0 = Tswitch(2);
        if t <= t0
            Kon = K0;
        else
            Kon = K0*exp(-m*(t-t0));
        end
end

%%%% This function generate ONE mysion contractile profile
    function [Fm] = OneMyosin(t, myo_alpha, f0, P)
        tloc = rem(t, P);    % Time modulo P (periodic behavior)
        a1 = myo_alpha(1);         % alpha_1
        a2 = myo_alpha(2);         % alpha_2
        a3 = myo_alpha(3);         % alpha_3

        % Creating the piecewise function for the myosin profile
        if (tloc <= a1*P)
            % Increasing phase
            Fm = f0 / (a1 * P) * tloc;
        elseif (tloc > a1*P) && (tloc <= a2*P)
            % Plateau phase
            Fm = f0;
        elseif (tloc > a2*P) && (tloc <= a3*P)
            % Decreasing phase
            Fm = -f0 / ((a3 - a2) * P) * (tloc - a2 * P) + f0;
        elseif (tloc > a3*P) && (tloc <= P)
            % Zero phase
            Fm = 0;
        end
    end

    
end