%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2024                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Approximation of velocity term in nonlocal 1D sensing region    %%%
%%%             [copyright: Chiara Villa (*), Nadia Loy (**)]           %%%
%%%                                                                     %%%
%%% (*) chiara.villa.1[at]sorbonne-universite.fr                        %%%
%%% (**) nadia.loy@polito.it                                            %%%
%%%                                                                     %%%
%%%   Disclaimer: this is not optimised (i.e. overall code will be slow %%%
%%%   if this needs to be re-calculated at every iteration in time),    %%%
%%%   for improvements see Gerisch 2010 (doi:10.1093/imanum/drp027)     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [UTx2A,UTx2B,DTx2A,DTx2B] = Nonlocal_advection_1D(par,Kappa)
%%%
%%% UT given at cell edges (not cell centers)
%%%
%%% For now: code given chosing gamma = Dirac delta cenetered in R(y)
%%% With options: Kappa either Dirac delta or Von mises distributions

a = par.x1min;
b = par.x1max;
c = par.x2min;
d = par.x2max;

% x1 and x2 discretisation at cell edges (maintaining same dx)
x1A = 0.5*(par.str1+a); % Middle of LM stripe
x1B = 0.5*(par.str2+par.str1); %Middle of FN stripe
x2 = (par.x2min):par.dx2:(par.x2max); 
Nx1 = 1; % 1D
Nx2 = length(x2);

% y discretisation (choose dy)
dy = 0.1;
y = (0:dy:1);
Ny = length(y);

% vhat discretisation
vhat = [-1,1];
th = [pi/2,3*pi/2];
dth = 1;

% ECM density gradient (determining direction of motion) - def (64)
Mat=@(s2) par.Mmin+par.Mgr*(s2-c);    
% Constant of integration (only if gamma = Dirac delta) > x2-dependent
cM = repmat(1./((2)*(par.Mmin+par.Mgr*(x2-c))),Nx1,1);  

% Mean post-reorientation speed - i.e. u_psi in def (57)
vbar_f=@(s1,s2) par.vmaxL*((s1<par.str1) | (s1 >par.str2)) + par.vmaxF*((s1>=par.str1) & (s1 <=par.str2));

% Mean post-transition phenotype - def (59)
% For the purpose of the 1D test: choose ybarLB for both and just evaluate
% A at LN stripe and B at FN stripe
ybar_fA=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));
ybar_fB=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));

% Sensing radius - def (55)
R_f=@(yd) par.Rmin + yd.*(par.Rmax-par.Rmin); 

% Initialise UT at cell edges
UTx2A = zeros(Nx1,Nx2);
UTx2B = zeros(Nx1,Nx2);
DTx2A = zeros(Nx1,Nx2);
DTx2B = zeros(Nx1,Nx2);

switch Kappa

    case "DD" % Dirac delta - Definition (60)

    % Evaluate ybar outside the loop if... ? (ybar determined locally?)
    ybA = repmat(ybar_fA(x1A,x2),1,Nx2);
    ybB = repmat(ybar_fB(x1B,x2),1,Nx2);
    
    for j=1:2
    
        % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
        % > avoid exiting the domain
        vbA=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1A+R_f(ybA).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
        vbB=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1B+R_f(ybB).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
        
        % Add integration contribution to UT
        UTx2A=UTx2A+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j)))))*sin(th(j)).*vbA*dth;
        UTx2B=UTx2B+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j)))))*sin(th(j)).*vbB*dth;
    
    end
    
    for j=1:2
    
        % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
        % > avoid exiting the domain
        vbA=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1A+R_f(ybA).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
        vbB=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1B+R_f(ybB).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
            
        % Add integration contribution to D_T 
        DTx2A=DTx2A+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))).*((vbA.*sin(th(j))-UTx2B).^2)*dth;
        DTx2B=DTx2B+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))).*((vbB.*sin(th(j))-UTx2B).^2)*dth;
    
    end

    case "VM" % Von Mises - Definition (61)

    k = par.ky; % Von Mises coefficient
    
    % Initialise also uT (required if Kappa is not a Dirac delta distribution)
    uTx2A=zeros(Nx1,Nx2);
    uTx2B=zeros(Nx1,Nx2);
    aT22A = zeros(Nx1,Nx2);
    DT22A = zeros(Nx1,Nx2);
    aT22B = zeros(Nx1,Nx2);
    DT22B = zeros(Nx1,Nx2);
       
    
    % Evaluate ybar outside the loop if... ?
    ybA=repmat(ybar_fA(x1',x2),1,Nx2);
    ybB=repmat(ybar_fB(x1',x2),1,Nx2);
    
    for i=1:Ny
        for j=1:Nth
    
            % Evaluate vbA at each (x1,x2) + R(y)*(cos(th),sin(th))
            % > avoid exiting the domain
            vbA = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+par.Rmax.*y(i).*cos(th(j)))),...
                min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+par.Rmax.*y(i).*sin(th(j)))));
            vbB = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+par.Rmax.*y(i).*cos(th(j)))),...
                min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+par.Rmax.*y(i).*sin(th(j)))));
    
            % Add integration contribution to uT
            uTx2A = uTx2A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                x2+par.Rmax.*y(i).*sin(th(j)))))*sin(th(j)).*vbA*dth;
            uTx2B = uTx2B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                x2+Rmax.*y(i).*sin(th(j)))))*sin(th(j)).*vbB*dth;

            % Add integration contribution to aT
            aT22A = aT22A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                x2+R_f(y(i)).*sin(th(j)))))*(sin(th(j)).^2).*(vbA.^2)*dth;
            aT22B = aT22B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                x2+R_f(y(i)).*sin(th(j)))))*(sin(th(j)).^2).*(vbB.^2)*dth;

        end
    
        % Evaluate UT integrating uT against Von Mises distribution
        UTx2A=UTx2A+uTx2A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)+2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
        UTx2B=UTx2B+uTx2B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)+2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;
 
        % Evaluate AT integrating aT against Von Mises distribution
        DT22A = DT22A + aT22A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
            + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
        DT22B = DT22B + aT22B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
            + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;

    end

    % Compute variance-covariace matrix DT = AT - UT⊗UT
    DTx2A = DT22A - UTx2A.^2;
    DTx2B = DT22A - UTx2A.^2;

otherwise

        error('Unadmissible definition of Kappa: either Dirac delta (DD) or Von Mises (VM)')
end


end               
