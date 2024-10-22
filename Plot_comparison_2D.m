%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Compare solutions to macroscopic PDE and mesoscopic kinetic model  %%%
%%%  of Stripe migration essay (1D cross section)                       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
folder = 'Saved_Data/';

%% Load kinetic model solution 
load([folder,'datasMC_2D_eps1e-3_N1e6.mat'])
Dx1=0.25;
Dx2=Dx1;
rho_micro_2D % function to reconstruct cell density

%% Load PDE solution
load([folder,'Saved_110924_Setup_DD_eps0.001.mat'])
load([folder,'Saved_110924_Results_DD_eps0.001.mat'])

Nx1 = (par.x1max-par.x1min)/par.dx1;
Nx2 = (par.x2max-par.x2min)/par.dx2;
x1 = (par.x1min+0.5*par.dx1):par.dx1:(par.x1max-0.5*par.dx1); 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); 

%%% Choose times at which to plot
T1=24;
T2=48;

rho24 = rhostoreB(:,T1*Nx2+1:(T1+1)*Nx2);
rho48 = rhostoreB(:,T2*Nx2+1:(T2+1)*Nx2);

%% Select cross-sections to plot longitudinally 

%%% Middle of first laminin stripe
x1plot_L = (par.x1min + par.str1)/2;
[~, index_L] = min(abs(x1 - x1plot_L));
x1plot_L = x1(index_L);

x1plot_Lc = a+(b-a)/6;
[~, index_Lc] = min(abs(xg1 - x1plot_Lc));
x1plot_Lc = xg1(index_Lc);

%%% Middle of central fibronectin stripe
x1plot_F = (par.str1 + par.str2)/2;
[~, index_F] = min(abs(x1 - x1plot_F));
x1plot_F = x1(index_F);

x1plot_Fc = a+(b-a)/2;
[~, index_Fc] = min(abs(xg1 - x1plot_Fc));
x1plot_Fc = xg1(index_Fc);

%% Plot
figure(1)

subplot(2,5,[1,2,3])
pbaspect([6 1 1])
xlim([par.x2min,par.x2max])
plot(x2,rho24(index_L,:),'Color', [0 0 1 0.5])
hold on
plot(x2,rho48(index_L,:),'Color', [0 0 1 1])

st=1;
plot(xg2(1:st:end-1),rhoc24(index_Lc,1:st:end),'o','Color', [1 0.5 0.5 0.5])
hold on
plot(xg2(1:st:end-1),rhoc48(index_Lc,1:st:end),'o', 'Color',[1 0 0 0.9])
legend('macro $t=24$','macro $t=48$','micro $t=24$','micro $t=48$','Interpreter','latex')
xlabel('$x_2$')
title(['Solution comparison at $x_1=$ ',num2str(x1(index_L)),' (over Laminin)'])

subplot(2,5,[6,7,8])
pbaspect([6 1 1])
xlim([par.x2min,par.x2max])
plot(x2,rho24(index_F,:),'Color', [0 0 1 0.5])
hold on
plot(x2,rho48(index_F,:),'Color', [0 0 1 1])

plot(xg2(1:end-1),rhoc24(index_Fc,:),'o','Color', [1 0.5 0.5 0.5])
hold on
plot(xg2(1:end-1),rhoc48(index_Fc,:),'o','Color', [1 0 0 0.8])
legend('macro $t=24$','macro $t=48$','micro $t=24$','micro $t=48$','Interpreter','latex')
xlabel('$x_2$')
title(['Solution comparison at $x_1=$ ',num2str(x1(index_F)),' (over Fibronectin)'],'Interpreter','latex')


%% Select cross-sections to plot transversally
x2p=-0.225;

x2plot1 = x2p;
[~, index_1] = min(abs(x2 - x2plot1));
x2plot1 = x2(index_1);

 x2plot_1c = x2p; %a+(b-a)/6;
 [~, index_1c] = min(abs(xg2 - x2plot_1c));
 x2plot_1c = xg2(index_1c);

%%% Plot
subplot(2,5,[4,5,9,10])
xlim([par.x1min,par.x1max])
plot(x1,rho24(:,index_1),'Color', [0 0 1 0.5])
hold on
plot(x1,rho48(:,index_1),'Color', [0 0 1 1])
legend(['PDE $t=$',num2str(T1)],['PDE $t=$',num2str(T2)],'Interpreter','latex')

%%% with kinetic solution 
st=1;
plot(xg1(1:st:end-1),rhoc24(1:st:end,index_1c),'o','Color', [1 0.5 0.5 0.5])
hold on
plot(xg1(1:st:end-1),rhoc48(1:st:end,index_1c),'o', 'Color',[1 0 0 0.9])
legend(['macro $t=$',num2str(T1)],['macro $t=$',num2str(T2)],['micro $t=$',num2str(T1)],['micro $t=$',num2str(T2)],'Interpreter','latex')

xlabel('$x_1$')
title(['Solution comparison at $x_2=$ ',num2str(x2(index_1))])

