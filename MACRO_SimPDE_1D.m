%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Modelling collective migration of phenotypically heterogeneous   %%%
%%%          cell populations: from single-cell dynamics                %%%     
%%%                to population-level behaviours"                      %%%
%%%                                                                     %%%
%%%            T. Lorenzi, N. Loy (*), C. Villa, 2024                   %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code to simulate Stripe migration essays of Goodman et al. (1989)  %%%
%%%  over laminin & fibronectin, with the macroscopic model of Eq.(46)  %%%
%%%                  [copyright: Chiara Villa (*)]                      %%%
%%%                                                                     %%%
%%% (*) chiara.villa.1[at]sorbonne-universite.fr                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
%%% To save data
date = '110924'; 
folder = 'Saved_Data/'; 
add = ''; % default

%% Parameters 

%%% Choose Kappa definition:
Kappa = 'DD'; % Dirac delta
% Kappa = 'VM'; % Von Mises

%%% Problem parameters (in units of 10^-2 cm)
par = Parameters();
par.eps = 0.001;
D = num2str(par.eps);

%%% Space discretisation
dx = 0.01;                   % Choose grid refinement
Nx1 = 1;
par.dx2 = dx; 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); % Cell centres
Nx2 = length(x2);

%%% Time discretisation
dt = 0.01; % NOTE: this will be updated by CFL condition after computing UT, DT
Tend = 48;
t = 0:dt:Tend;
Nt = length(t);
tstore = [0:1:Tend];

%% Initial conditions
rho0 = zeros(1,Nx2);
L2 = 1;
rho0(1,x2<-L2) = 1;
rho0(1,x2>=-L2) = 1-(((x2(x2>=-L2)+L2)./L2).^2);
rho0(1,x2<2*L2-par.x2max) = 1-(((2*L2-par.x2max-x2(x2<2*L2-par.x2max))./L2).^2);
rho0(rho0<0) = 0;
% Calculate initial mass to check mass conservation
Mass = sum(rho0,"all")*par.dx2; 

%% Computation of the advection velocity
%  > Details of the problem-specific definitions affecting UT are given in 
%    the file 'Nonlocal_advection.m' 
%  > UT and DT are calculated and given at the cell edges

% 1D: A = LN stripe and B = FN stripe
[UTx2A,UTx2B,DTx2A,DTx2B] = Nonlocal_advection_1D(par,Kappa);

save([folder,'Saved_',date,'_Setup_1D_eps',D,'.mat'],'par','Kappa','UTx2A',"UTx2B",'DTx2A',"DTx2B")

%% Setup solver for equation (70)-(71) - Advection velocity UT(1-eps\div UT)
%%%
%%% \dt p + \div ( p * UT(1-eps\div UT) ) = \div \div (eps * DT * p )

%%% Ghost points to compute \div UT (ensuring no-flux boundary conditions)
UTx2gpA = [UTx2A(:,1),UTx2A,UTx2A(:,end)];
UTx2gpB = [UTx2B(:,1),UTx2B,UTx2B(:,end)];

%%% First order correction contribution to advection velocity
UTx2A = UTx2A.*( 1 - par.eps*((UTx2gpA(:,3:end)-UTx2gpA(:,1:end-2))./(2*par.dx2)) );
UTx2B = UTx2B.*( 1 - par.eps*((UTx2gpB(:,3:end)-UTx2gpB(:,1:end-2))./(2*par.dx2)) );

%%% D_T needs to be given at cell interface
DTA_c = 0.5*(DTx2A(2:end)+DTx2A(1:end-1));
DTB_c = 0.5*(DTx2B(2:end)+DTx2B(1:end-1));

%%% Explicit diffusion matrix: eps * \div \div (...)
DMx2 = DM_def(Nx2,par.dx2,par.eps);

%% Iterate in time (MUSCL scheme)

%%% Store 
rhoA = rho0;
rhostoreA = [rhoA];
rhoB = rho0;
rhostoreB = [rhoB];

%%% CFL condition
UmaxA = max(max(max(UTx1A)),max(max(UTx2A)));
UmaxB = max(max(max(UTx1B)),max(max(UTx2B)));
Umax = max(UmaxA,UmaxB);
DmaxA = max(max(max(DTx1A)),max(max(DTx2A)));
DmaxB = max(max(max(DTx1B)),max(max(DTx2B)));
Dmax = par.eps*max(DmaxA,DmaxB);
CFL_adv = dx/Umax;
CFL_diff = (dx^2)/(2*Dmax);
if dt>CFL_adv || dt>CFL_diff
   dt = 0.5*min(CFL_adv,CFL_diff);
   % CFL condition changed the dt, ensure solutions will be stored hourly
   for i=floor(1.0/dt):-1:1
       if mod(1,1.0/i)==0
           dt = 1.0/i
           break;
       end
   end
end

%%% Iterate

for i=dt:dt:Tend

    %%% Compute df/dx given f=UT*rho, using MUSCL: add ghost points
    rhogp2A = [zeros(size(rhoA,1),2),rhoA,zeros(size(rhoA,1),2)];
    dfdx2A = MUSCL_GP(rhogp2A,UTx2A,Nx2,3,par.dx2,dt);

    rhogp2B = [zeros(size(rhoB,1),2),rhoB,zeros(size(rhoB,1),2)];
    dfdx2B = MUSCL_GP(rhogp2B,UTx2B,Nx2,3,par.dx2,dt);

    %%% Explicit in time
    rhoA = rhoA - dt*(dfdx2A) + dt*(DMx2*(DTA_c.*(rhoA))')';

    rhoB = rhoB - dt*(dfdx2B) + dt*(DMx2*(DTB_c.*(rhoB))')';

    %%% Store Mass to check Mass Conserrvation
    Mass = [Mass,sum(rhoB,"all")*par.dx2];

    %%% Store every hour
    if ismember(round(i,6),tstore)%mod(1.0,i)==0
        rhostoreA = [rhostoreA,rhoA];
        rhostoreB = [rhostoreB,rhoB];
        subplot(2,3,1)
        plot(UTtest)
        title('UT original')
        subplot(2,3,2)
        plot(UTx2A)
        title('UTe with eps correction')
        subplot(2,3,3)
        plot(dfdx2A )
        title('advection d(UTe*rho)/dx')
        subplot(2,3,4)
        plot(DTA_c)
        title('DT')
        subplot(2,3,5)
        plot(d2fdx2A )
        title('diffusion d2(DT*rho)/dx2')
        subplot(2,3,6)
        plot(rhoA )
        title('rho')
        drawnow
     end

end

save([folder,'Saved_',date,'_Results_1D_eps',D,'.mat'],'rhostoreA','rhostoreB','Mass','rhoA','rhoB','rho0')

%% Final plot
% load(['Saved_',date,'_Setup.mat'])
% load(['Saved_',date,'_Results.mat'])

%%% Plot stripes
figure(1)
subplot(1,2,1)
plot(x2,rho0,'k:')
hold on
plot(x2,rhoA,'r')
title('LN')
subplot(1,2,2)
plot(x2,rho0,'k:')
hold on
plot(x2,rhoB,'b')
title('FN')

%%% Plot to check mass conservation
figure(2)
plot(Mass)

%% Diffusion matrix (finite volume scheme with no-flux BCs)

function DM = DM_def(Nx,dx,D)

    %%% Coefficient matrix
    DM = -2*eye(Nx);
    DM(2:Nx,1:Nx-1) = DM(2:Nx,1:Nx-1) + eye(Nx-1);
    DM(1:Nx-1,2:Nx) = DM(1:Nx-1,2:Nx) + eye(Nx-1);

    %%% No-flux boundary conditions (from finite volume scheme)
    DM(1,1) = -1;
    DM(1,2) = 1;
    DM(Nx,Nx-1) = 1;
    DM(Nx,Nx) = -1;

    %%% Output
    DM = D*DM./(dx^2);

end
