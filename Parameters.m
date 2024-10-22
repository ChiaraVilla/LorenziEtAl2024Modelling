%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Input parameters for simulated problem (not numerical parameters)  %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par = Parameters() % In units of 10^(-2) cm, i.e. 100 micrometers

par.vmaxL = 4*0.5;   % Maximum laminin-induced post-reorientation speed
par.vmaxF = 4*0.15;  % Maximum fibronectin-induced post-reorientation speed
par.ymaxF = 0;       % Maximum fibronectin-induced post-transition phenotype

% Maximum laminin-induced post-transition phenotype
par.ymaxLA = 0;      % Option A
par.ymaxLB = 1;      % Option B

par.sig2 = 0.01;     % Variance of phenotypic transition probability kernel
par.Rmax = 0.5;      % Maximum cell sensing radius (0.005 cm)
par.Rmin = 0.05;     % Minimum cell sensing radius (0.0005 cm)
par.L0 = 1;          % Baseline constant laminin                                  
par.F0 = 1;          % Baseline constant fibronectin    
par.eps = 1e-4;      % First order correction coefficient 
par.ky = 1;          % Von Mises coefficient

% Domain bounds
par.x1min = 0;         
par.x1max = 6;
par.x2min = -6;
par.x2max = 6;

% Stripe boundaries (3 stripes)
par.str1 = (par.x1max-par.x1min)/3.0;
par.str2 = 2*(par.x1max-par.x1min)/3.0;

% ECM density parameters
par.Mmin = 0.1;
par.Mgr = 1;

end