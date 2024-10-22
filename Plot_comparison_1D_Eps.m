%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Compare solutions to macroscopic PDE and mesoscopic micro model  %%%
%%%  of Stripe migration essay (1D cross section)                       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
folder = 'Saved_Data/';

%% Load PDE solution
epsval=10^-4;
load([folder,'datasMC_1D_tot_eps',num2str(epsval),'.mat'])
rho_micro1D
load([folder,'Saved_11092024_Results_1D_eps',num2str(epsval),'.mat'])
load([folder,'Saved_11092024_Setup_1D_eps',num2str(epsval),'.mat'])

Nx1 = 1;
Nx2 = (par.x2max-par.x2min)/par.dx2;
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); 


%% Plot
% CHOOSE AT WHICH TIMES TO PLOT
T1 = 12;
T2 = 24;

rho12 = rhostoreA(:,T1*Nx2+1:(T1+1)*Nx2);
rho24 = rhostoreA(:,T2*Nx2+1:(T2+1)*Nx2);

subplot(4,1,1)
pbaspect([6 1 1])
xlim([par.x2min,par.x2max])
%%% PDE solution
plot(x2,rho0,'Color', [0 0 1 0.2])
hold on
plot(x2,rho12,'Color', [0 0 1 0.5])
hold on
plot(x2,rho24,'Color', [0 0 1 1])
%%% micro solution 
st=1;
plot(xg2(1:st:end-1),rhoc0(1:st:end),'o','Color', [1 0.7 0.7 0.2])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc12(1:st:end),'o','Color', [1 0.5 0.5 0.5])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc24(1:st:end),'o', 'Color',[1 0 0 1])%,'HandleVisibility','off')
legend('macro $t=0$',['macro $t=$',num2str(T1)],['macro $t=$',num2str(T2)],'micro $t=0$',['micro $t=$',num2str(T1)],['micro $t=$',num2str(T2)],'Interpreter','latex')
xlabel('$x_2$')
title('$\varepsilon=10^{-4}$')
ylim([0 3.2])


%% Load PDE solution
epsval = 1e-3;
load([folder,'datasMC_1D_tot_eps',num2str(epsval),'.mat'])
rho_micro1D
load([folder,'Saved_11092024_Results_1D_eps',num2str(epsval),'.mat'])
load([folder,'Saved_11092024_Setup_1D_eps',num2str(epsval),'.mat'])

Nx1 = 1;
Nx2 = (par.x2max-par.x2min)/par.dx2; 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); 

%% Plot
% CHOOSE AT WHICH TIMES TO PLOT
T1 = 12;
T2 = 24;

rho12 = rhostoreA(:,T1*Nx2+1:(T1+1)*Nx2);
rho24 = rhostoreA(:,T2*Nx2+1:(T2+1)*Nx2);

subplot(4,1,2)
pbaspect([6 1 1])
xlim([par.x2min,par.x2max])
%%% PDE solution
plot(x2,rho0,'Color', [0 0 1 0.2])
hold on
plot(x2,rho12,'Color', [0 0 1 0.5])
hold on
plot(x2,rho24,'Color', [0 0 1 1])
%%% micro solution 
st=1;
plot(xg2(1:st:end-1),rhoc0(1:st:end),'o','Color', [1 0.7 0.7 0.2])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc12(1:st:end),'o','Color', [1 0.5 0.5 0.5])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc24(1:st:end),'o', 'Color',[1 0 0 1])%,'HandleVisibility','off')
legend('macro $t=0$',['macro $t=$',num2str(T1)],['macro $t=$',num2str(T2)],'micro $t=0$',['micro $t=$',num2str(T1)],['micro $t=$',num2str(T2)],'Interpreter','latex')
xlabel('$x_2$')
title('$\varepsilon=10^{-3}$')
ylim([0 3.2])

%% Load PDE solution
epsval = 1e-2;
load([folder,'datasMC_1D_tot_eps',num2str(epsval),'.mat'])
rho_micro1D
load([folder,'/Saved_11092024_Results_1D_eps',num2str(epsval),'.mat'])
load([folder,'Saved_11092024_Setup_1D_eps',num2str(epsval),'.mat'])

Nx1 = 1;
Nx2 = (par.x2max-par.x2min)/par.dx2;
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); 

%% Plot
% CHOOSE AT WHICH TIMES TO PLOT
T1 = 12;
T2 = 24;

rho12 = rhostoreA(:,T1*Nx2+1:(T1+1)*Nx2);
rho24 = rhostoreA(:,T2*Nx2+1:(T2+1)*Nx2);

subplot(4,1,3)
pbaspect([6 1 1])
xlim([par.x2min,par.x2max])
%%% PDE solution
plot(x2,rho0,'Color', [0 0 1 0.2])
hold on
plot(x2,rho12,'Color', [0 0 1 0.5])
hold on
plot(x2,rho24,'Color', [0 0 1 1])
%%% micro solution 
st=1;
plot(xg2(1:st:end-1),rhoc0(1:st:end),'o','Color', [1 0.7 0.7 0.2])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc12(1:st:end),'o','Color', [1 0.5 0.5 0.5])%,'HandleVisibility','off')
hold on
plot(xg2(1:st:end-1),rhoc24(1:st:end),'o', 'Color',[1 0 0 1])%,'HandleVisibility','off')
legend('macro $t=0$',['macro $t=$',num2str(T1)],['macro $t=$',num2str(T2)],'micro $t=0$',['micro $t=$',num2str(T1)],['micro $t=$',num2str(T2)],'Interpreter','latex')
xlabel('$x_2$')
title('$\varepsilon=10^{-2}$')
ylim([0 3.2])
