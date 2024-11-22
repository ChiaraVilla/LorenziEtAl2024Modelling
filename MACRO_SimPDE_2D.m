%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2024                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code to simulate Stripe migration essays of Goodman et al. (1989)  %%%
%%%  over laminin & fibronectin, with the macroscopic model of Eq.(44)  %%%
%%%                 in 2D [copyright: Chiara Villa (*)]                 %%%
%%%                                                                     %%%
%%% (*) chiara.villa.1[at]sorbonne-universite.fr                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 10)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
%%% To save data
date = '110924'; 
folder = 'Saved_Data/'; 
add = ''; % default

%% Parameters and discretisation

%%% Choose Kappa definition:
% Kappa = 'DD'; % Dirac delta
Kappa = 'VM'; % Von Mises

%%% Problem parameters (in units of 10^-2 cm)
par = Parameters();
par.eps = 0.001;
add = ['_eps',num2str(par.eps)]; % any additional specification for data file titles

%%% Space discretisation
dx = 0.05; % Choose grid refinement
par.dx1 = dx; 
x1 = (par.x1min+0.5*par.dx1):par.dx1:(par.x1max-0.5*par.dx1); % Cell centres
Nx1 = length(x1);
par.dx2 = dx; 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); % Cell centres
Nx2 = length(x2);

%%% Final time
Tend = 48;

%% Initial conditions
rho0 = zeros(Nx1,Nx2);
L2 = 1;
rho0(1:Nx1,x2<-L2) = 1;
for i=1:Nx1
    rho0(i,x2>=-L2) = 1-(((x2(x2>=-L2)+L2)./L2).^2);
    rho0(i,x2<2*L2-par.x1max) = 1-(((2*L2-par.x1max-x2(x2<2*L2-par.x1max))./L2).^2);
end
rho0(rho0<0) = 0;
% Calculate initial mass to check mass conservation
MassA = sum(rho0,"all")*par.dx2*par.dx1; 
MassB = MassA;

%% Computation of the advection velocity
%  > Details of the problem-specific definitions affecting UT are given in 
%    the file 'Nonlocal_advection.m' 
%  > UT is calculated and given at the cell edges

[UTx1A,UTx2A,UTx1B,UTx2B,DT11A,DT12A,DT22A,DT11B,DT12B,DT22B] = Nonlocal_advection_2D(par,Kappa);

save([folder,'Saved_',date,'_Setup_',Kappa,add,'.mat'],'par','Kappa','UTx2A',"UTx1A","UTx1B","UTx2B",'DT11A','DT12A','DT22A','DT11B','DT12B','DT22B')

%% Setup solver for equation (73)-(74) - Advection velocity UT(1-eps\div UT)
%%%
%%% \dt p + \div ( p * UT(1-eps\div UT) ) = \div \div (eps * DT * p )

%%% Ghost points to compute \div UT (ensuring no-flux boundary conditions)
UTx1gpA = [UTx1A(1,:);UTx1A;UTx1A(end,:)];
UTx2gpA = [UTx2A(:,1),UTx2A,UTx2A(:,end)];
UTx1gpB = [UTx1B(1,:);UTx1B;UTx1B(end,:)];
UTx2gpB = [UTx2B(:,1),UTx2B,UTx2B(:,end)];

%%% First order correction contribution to advection velocity
UTx1A = UTx1A.*( 1 - par.eps*((UTx1gpA(3:end,:)-UTx1gpA(1:end-2,:))./(2*par.dx1)...
    + (UTx2gpA(:,3:end)-UTx2gpA(:,1:end-2))./(2*par.dx2)) );
UTx2A = UTx2A.*( 1 - par.eps*((UTx1gpB(3:end,:)-UTx1gpB(1:end-2,:))./(2*par.dx1)...
    + (UTx2gpB(:,3:end)-UTx2gpB(:,1:end-2))./(2*par.dx2)) );

%%% Explicit diffusion matrix
DMx1 = DM_def(Nx1,par.dx1,par.eps);
DMx2 = DM_def(Nx2,par.dx2,par.eps);

%% Iterate in time (MUSCL scheme)

%%% Store 
rhoA = rho0;
rhostoreA = [rhoA];
rhoB = rho0;
rhostoreB = [rhoB];

%%% Time discretisation & the CFL condition
dt = 0.01; % Choose maximum dt
UmaxA = max(max(max(UTx1A)),max(max(UTx2A)));
UmaxB = max(max(max(UTx1B)),max(max(UTx2B)));
Umax = max(UmaxA,UmaxB);
DmaxA = max(max(max(max(DT11A)),max(max(DT22A))),max(max(DT12A)));
DmaxB = max(max(max(max(DT11B)),max(max(DT22A))),max(max(DT12B)));
Dmax = max(DmaxA,DmaxB);
CFL_adv = dx/Umax;
CFL_diff = (dx^2)/(2*Dmax);
if dt>CFL_adv || dt>CFL_diff
   dt = 0.5*min(CFL_adv,CFL_diff);
   % CFL condition changed the dt, ensure solutions will be stored hourly
   for i=floor(1.0/dt):-1:1
       if mod(1,1.0/i)==0
           dt = 1.0/i;
           break;
       end
   end
end
t = 0:dt:Tend;
Nt = length(t);

%%% Time iteration
for i=1:Nt-1

    %%% Compute df/dx given f=UT*rho, using MUSCL: add ghost points
    rhogp1A = [zeros(2,size(rhoA,2));rhoA;zeros(2,size(rhoA,2))];
    dfdx1A = MUSCL_GP(rhogp1A',UTx1A',Nx1,3,par.dx1,dt)';
    rhogp2A = [zeros(size(rhoA,1),2),rhoA,zeros(size(rhoA,1),2)];
    dfdx2A = MUSCL_GP(rhogp2A,UTx2A,Nx2,3,par.dx2,dt);

    rhogp1B = [zeros(2,size(rhoB,2));rhoB;zeros(2,size(rhoB,2))];
    dfdx1B = MUSCL_GP(rhogp1B',UTx1B',Nx1,3,par.dx1,dt)';
    rhogp2B = [zeros(size(rhoB,1),2),rhoB,zeros(size(rhoB,1),2)];
    dfdx2B = MUSCL_GP(rhogp2B,UTx2B,Nx2,3,par.dx2,dt);

    %%% PDE
    rhoA = rhoA - dt*(dfdx1A+dfdx2A) + dt*( DMx1*(DT11A.*rhoA) + ...
        (DMx2*(DT22A.*rhoA)')' + 2*Ddiag(Nx1,Nx2,par.dx1,par.dx2,par.eps,DT12A.*rhoA) );

    rhoB = rhoB - dt*(dfdx1B+dfdx2B) + dt*( DMx1*(DT11B.*rhoB) + ...
        (DMx2*(DT22B.*rhoB)')' + 2*Ddiag(Nx1,Nx2,par.dx1,par.dx2,par.eps,DT12B.*rhoB) );

    %%% Store Mass to check Mass Conserrvation
    MassA = [MassA,sum(rhoA,"all")*par.dx2*par.dx1];
    MassB = [MassB,sum(rhoB,"all")*par.dx2*par.dx1];

    %%% Store every hour
    if mod(i-1,1/dt)==0
        rhostoreA = [rhostoreA,rhoA];
        rhostoreB = [rhostoreB,rhoB];
        subplot(1,2,1)
        Plot_StripeS_Timet(rhoA,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,max(max(rhoA)),'n')
        title(['$t=$',num2str(t(i))],'Interpreter','latex')
        subplot(1,2,2)
        Plot_StripeS_Timet(rhoB,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,max(max(rhoB)),'n')
        title(['$t=$',num2str(t(i))],'Interpreter','latex')
        drawnow
    end

end

save([folder,'Saved_',date,'_Results_',Kappa,add,'.mat'],'rhostoreA','rhostoreB','MassA','MassB','rhoA','rhoB','rho0')

%% Final plot
% Kappa = 'DD';
% load([folder,'Saved_',date,'_Setup_',Kappa,add,'.mat'])
% load([folder,'Saved_',date,'_Results_',Kappa,add,'.mat'])

%%% Plot Stripes
figure(1)
Plot_stripes(rho0,rhoA,rhoB,par,x1,x2)

%%% Plot to check mass conservation
figure(2)
subplot(1,2,1)
plot(MassA)
title('Mass $\rho_A$','Interpreter','latex')
subplot(1,2,2)
plot(MassB)
title('Mass $\rho_B$','Interpreter','latex')

%% Diffusion matrix (finite volume scheme with no-flux BCs)

function DM = DM_def(Nx,dx,D) % 1D Fickian diffusion (d_xx) 

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

function Dd = Ddiag(Nx1,Nx2,dx1,dx2,D,A) % Second order mixed derivative

    %%% Initialize matrix 
    Dd = zeros(Nx1,Nx2);
    
    %%% Internal entries 
    Dd(2:Nx1-1,2:Nx2-1) = ( A(3:Nx1,3:Nx2) - A(3:Nx1,1:Nx2-2) ...
                           - A(1:Nx1-2,3:Nx2) + A(1:Nx1-2,1:Nx2-2) );
    
    %%% Ghost points to ensure no-flux BCs are satisfied
    A_l = A(1,:);     % left       
    A_r = A(Nx1,:);   % right 
    A_b = A(:,1);     % bottom      
    A_t = A(:,Nx2);   % top    
    A_lb = A(1,1);
    A_rb = A(Nx1,1);
    A_lt = A(1,Nx2);
    A_rt = A(Nx1,Nx2);
    
    %%% Boundaries
    Dd(1,2:Nx2-1) = A(2,3:Nx2) - A(2,1:Nx2-2) - A_l(3:Nx2) + A_l(1:Nx2-2);
    Dd(Nx1,2:Nx2-1) = A_r(3:Nx2) - A_r(1:Nx2-2) - A(Nx1-1,3:Nx2) + A(Nx1-1,1:Nx2-2);
    Dd(2:Nx1-1,1) = A(3:Nx1,2) - A(1:Nx1-2,2) - A_b(2:Nx1-1) + A_b(1:Nx1-2);
    Dd(2:Nx1-1,Nx2) = A(3:Nx1,Nx2-1) - A(1:Nx1-2,Nx2-1) - A_t(2:Nx1-1) + A_t(1:Nx1-2);
    Dd(1,1) = A(2,2) - A_b(2) - A_l(2) + A_lb;
    Dd(1,Nx2) = A_t(2) - A(2,Nx2-1) + A_l(Nx2-1) - A_lt;
    Dd(Nx1,1) = A_r(2) - A_rb + A_b(Nx1-1) - A(Nx1-1,2);
    Dd(Nx1,Nx2) = A_rt - A_r(Nx2-1) - A_t(Nx1-1) + A(Nx1-1,Nx2-1);
    
    %%% Output
    Dd = D*Dd./(4*dx1*dx2);

end

