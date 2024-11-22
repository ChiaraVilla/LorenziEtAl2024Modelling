%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2024                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Approximation of velocity term in nonlocal 2D sensing region    %%%
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

function [UTx1A,UTx2A,UTx1B,UTx2B,DT11A,DT12A,DT22A,DT11B,DT12B,DT22B] = Nonlocal_advection_2D(par,Kappa)
%%%
%%% UT given at cell edges, DT given at cell centers
%%%
%%% Code given chosing: 
%%% - gamma(r) = Dirac delta centered in R(y) > definition (10)
%%% - psi(v) = Dirac delta centered in u_\psi > satisfies (11)
%%% - M(x) = Mmin + Mgr( x2 - Lm )
%%%
%%% Kappa(y): Dirac delta (DD) or Von Mises (VM) distribution

a = par.x1min;
b = par.x1max;
c = par.x2min;
d = par.x2max;

% x1 and x2 discretisation at cell edges (maintaining same dx)
x1 = (par.x1min):par.dx1:(par.x1max);
x2 = (par.x2min):par.dx2:(par.x2max); 
Nx1 = length(x1);
Nx2 = length(x2);

% y discretisation (choose dy)
dy = 0.1;
y = (0:dy:1);
Ny = length(y);

% theta discretisation (choose Nth to ensure equally spaced points)
Nth = 60;
dth = 2*pi/Nth;
th = (0:dth:2*pi-dth)';

% ECM density gradient (determining direction of motion) - def (64)
Mat=@(s2) par.Mmin+par.Mgr*(s2-c);    
% Constant of integration (only if gamma = Dirac delta) > x2-dependent   
cM = repmat(1./((2*pi)*(par.Mmin+par.Mgr*(x2-c))),Nx1,1);    

% Mean post-reorientation speed - i.e. u_psi in def (57)
vbar_f=@(s1,s2) par.vmaxL*((s1<par.str1) | (s1 >par.str2)) + par.vmaxF*((s1>=par.str1) & (s1 <=par.str2));

% Mean post-transition phenotype - def (59)
ybar_fA=@(s1,s2) par.ymaxLA*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));
ybar_fB=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));

% Sensing radius - def (55)
R_f=@(yd) par.Rmin + yd.*(par.Rmax-par.Rmin); 

% Initialise UT, AT and DT at cell edges
UTx1A = zeros(Nx1,Nx2);
UTx2A = zeros(Nx1,Nx2);
DT11A_i = zeros(Nx1,Nx2);
DT12A_i = zeros(Nx1,Nx2);
DT22A_i = zeros(Nx1,Nx2);

UTx1B = zeros(Nx1,Nx2);
UTx2B = zeros(Nx1,Nx2);
DT11B_i = zeros(Nx1,Nx2);
DT12B_i = zeros(Nx1,Nx2);
DT22B_i = zeros(Nx1,Nx2);

switch Kappa

    case "DD" % Dirac delta

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate ybar at each (x1,x2)
        ybA = repmat(ybar_fA(x1',x2),1,Nx2);
        ybB = repmat(ybar_fB(x1',x2),1,Nx2);
        
        for j=1:Nth
        
            % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
            % > avoid exiting the domain
            vbA = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(ybA).*cos(th(j)))),...
                min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
            vbAs(j,1:Nx1,1:Nx2) = vbA;
        
            vbB = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(ybB).*cos(th(j)))),...
                min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
            vbBs(j,1:Nx1,1:Nx2) = vbB;

            % Evaluate integrand contribution from ECM
            cMA = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j)))));
            cMAs(j,1:Nx1,1:Nx2) = cMA;

            cMB = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j)))));
            cMBs(j,1:Nx1,1:Nx2) = cMB;
        
            % Add integration contribution to UT
            UTx1A = UTx1A + cMA*cos(th(j)).*vbA*dth;
            UTx2A = UTx2A + cMA*sin(th(j)).*vbA*dth;
        
            UTx1B = UTx1B + cMB*cos(th(j)).*vbB*dth;
            UTx2B = UTx2B + cMB*sin(th(j)).*vbB*dth;
        
        end

        for j=1:Nth

            cMA(1:Nx1,1:Nx2) = cMAs(j,1:Nx1,1:Nx2);
            vbA(1:Nx1,1:Nx2) = vbAs(j,1:Nx1,1:Nx2);

            cMB(1:Nx1,1:Nx2) = cMBs(j,1:Nx1,1:Nx2);
            vbB(1:Nx1,1:Nx2) = vbBs(j,1:Nx1,1:Nx2);

            % Add integration contribution to DT
            DT11A_i = DT11A_i + cMA.*((vbA.*cos(th(j))-UTx1A).^2)*dth;
            DT12A_i = DT12A_i + cMA.*(vbA.*cos(th(j))-UTx1A).*(vbA.*sin(th(j))-UTx2A)*dth;
            DT22A_i = DT22A_i + cMA.*((vbA.*sin(th(j))-UTx2A).^2)*dth;
        
            DT11B_i = DT11B_i + cMB.*((vbB.*cos(th(j))-UTx1B).^2)*dth;
            DT12B_i = DT12B_i + cMB.*(vbB.*cos(th(j))-UTx1B).*(vbB.*sin(th(j))-UTx2B)*dth;
            DT22B_i = DT22B_i + cMB.*((vbB.*sin(th(j))-UTx2B).^2)*dth;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case "VM" % Von Mises

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        k = par.ky; % Von Mises coefficient
        
        % Initialise uT and aT (if Kappa is not a Dirac delta distribution)
        uTx1A = UTx1A;
        uTx2A = UTx2A;
        aT11A = zeros(Nx1,Nx2);
        aT12A = zeros(Nx1,Nx2);
        aT22A = zeros(Nx1,Nx2);
        DT11A = zeros(Nx1,Nx2);
        DT12A = zeros(Nx1,Nx2);
        DT22A = zeros(Nx1,Nx2);

        uTx1B = UTx1B;
        uTx2B = UTx2B;
        aT11B = zeros(Nx1,Nx2);
        aT12B = zeros(Nx1,Nx2);
        aT22B = zeros(Nx1,Nx2);
        DT11B = zeros(Nx1,Nx2);
        DT12B = zeros(Nx1,Nx2);
        DT22B = zeros(Nx1,Nx2);
        
        % Evaluate ybar at each (x1,x2)
        ybA=repmat(ybar_fA(x1',x2),1,Nx2);
        ybB=repmat(ybar_fB(x1',x2),1,Nx2);
        
        for i=1:Ny
            for j=1:Nth
        
                % Evaluate vbA at each (x1,x2) + R(y)*(cos(th),sin(th))
                % > avoid exiting the domain
                vbA = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(y(i)).*cos(th(j)))),...
                    min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(y(i)).*sin(th(j)))));
        
                vbB = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(y(i)).*cos(th(j)))),...
                    min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(y(i)).*sin(th(j)))));
        
                % Add integration contribution to uT
                uTx1A = uTx1A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*cos(th(j)).*vbA*dth;
                uTx2A = uTx2A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*sin(th(j)).*vbA*dth;
        
                uTx1B = uTx1B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*cos(th(j)).*vbB*dth;
                uTx2B = uTx2B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*sin(th(j)).*vbB*dth;

                % Add integration contribution to aT
                aT11A = aT11A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(cos(th(j)).^2).*(vbA.^2)*dth;
                aT12A = aT12A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(cos(th(j)).*sin(th(j))).*(vbA.^2)*dth;
                aT22A = aT22A + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(sin(th(j)).^2).*(vbA.^2)*dth;
        
                aT11B = aT11B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(cos(th(j)).^2).*(vbB.^2)*dth;
                aT12B = aT12B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(cos(th(j)).*sin(th(j))).*(vbB.^2)*dth;
                aT22B = aT22B + cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
                    x2+R_f(y(i)).*sin(th(j)))))*(sin(th(j)).^2).*(vbB.^2)*dth;
            end
        
            % Evaluate UT integrating uT against Von Mises distribution
            UTx1A = UTx1A + uTx1A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
                + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
            UTx2A = UTx2A + uTx2A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
                + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;

            UTx1B = UTx1B + uTx1B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
                + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;
            UTx2B = UTx2B + uTx2B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
                + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;

            % Evaluate AT integrating aT against Von Mises distribution
            DT11A = DT11A + aT11A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
                + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
            DT12A = DT12A + aT12A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
                + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
            DT22A = DT22A + aT22A.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
                + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;

            DT11B = DT11B + aT11B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
                + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;
            DT12B = DT12B + aT12B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
                + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;
            DT22B = DT22B + aT22B.*(0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
                + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;

        end

        % Compute variance-covariace matrix DT = AT - UT⊗UT
        DT11A_i = DT11A - UTx1A.^2;
        DT12A_i = DT12A - UTx1A.*UTx2A;
        DT22A_i = DT22A - UTx2A.^2;

        DT11B_i = DT11A - UTx1A.^2;
        DT12B_i = DT12A - UTx1A.*UTx2A;
        DT22B_i = DT22A - UTx2A.^2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    otherwise

        error('Unadmissible definition of Kappa: either Dirac delta (DD) or Von Mises (VM)')
end


% Provide DT at cell centers, not cell interfaces
DT11A = 0.25*(DT11A_i(2:end,2:end)+DT11A_i(1:end-1,2:end)+DT11A_i(2:end,1:end-1)+DT11A_i(1:end-1,1:end-1));
DT12A = 0.25*(DT12A_i(2:end,2:end)+DT12A_i(1:end-1,2:end)+DT12A_i(2:end,1:end-1)+DT12A_i(1:end-1,1:end-1));
DT22A = 0.25*(DT22A_i(2:end,2:end)+DT22A_i(1:end-1,2:end)+DT22A_i(2:end,1:end-1)+DT22A_i(1:end-1,1:end-1));

DT11B = 0.25*(DT11B_i(2:end,2:end)+DT11B_i(1:end-1,2:end)+DT11B_i(2:end,1:end-1)+DT11B_i(1:end-1,1:end-1));
DT12B = 0.25*(DT12B_i(2:end,2:end)+DT12B_i(1:end-1,2:end)+DT12B_i(2:end,1:end-1)+DT12B_i(1:end-1,1:end-1));
DT22B = 0.25*(DT22B_i(2:end,2:end)+DT22B_i(1:end-1,2:end)+DT22B_i(2:end,1:end-1)+DT22B_i(1:end-1,1:end-1));

end               
