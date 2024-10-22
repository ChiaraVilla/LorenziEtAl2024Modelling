# LorenziEtAl2024Modelling
Matlab code used to simulate the microscopic and macroscopic models in "Modelling collective migration of phenotypically heterogeneous cell populations: from single-cell dynamics to population-level behaviours" (2024), by Tommaso Lorenzi, Nadia Loy, Luigi Preziosi, Chiara Villa

## Generalities

**Public gitlab repository LorenziEtAl2024Modelling** <br />
This repository provides Matlab files to simulate the microscopic and macroscopic models in <br />
Tommaso Lorenzi, Nadia Loy, Luigi Preziosi, Chiara Villa (2024) <br />
Modelling collective migration of phenotypically heterogeneous cell populations: from single-cell dynamics to population-level behaviours <br />
Soon available on ArXiv [...] and HAL [...] <br />
For details we refer the interested reader to this publication. 

**Authors** <br />
Nadia Loy (Politecnico di Torino) and Chiara Villa (Sorbonne Université)

**Citation** <br />
Loy, N. and Villa, C. (2024). Matlab code to simulate microscopic and macroscopic models of collective cell migration on phenotypically heterogeneous cell populations. DOI soon available [...] <br />
If you use this software in your work then please cite the above named paper.

**Copyright notice** <br />
Matlab code to simulate microscopic and macroscopic models of collective cell migration on phenotypically heterogeneous cell populations. <br />
Copyright (C) 2024 N. Loy & C. Villa

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/.


## Repository content and how to use

**To change the model set up:** <br />
- 'Parameters.m' : the main parameter values of the simulations can be modified here, including the value of the first order correction coefficient 'epsilon'. This does not include numerical scheme-specific parameters, which can be modified in the respective files for simulating the macroscopic (SimPDE) and microscopic (SimMC) models. This function is called from the files:  MICRO_SimMC_1D.m, MICRO_SimMC_1D_tot.m, MICRO_SimMC_2D.m, MACRO_SimPDE_1D.m, MACRO_SimPDE_2D.m. 


### Simulating the microscopic model

The code is set up to simulate the microscopic model (1)-(2) described in Section 2 of LorenziEtAl2024Modelling, under the numerical setup detailed in Section 4.1 of LorenziEtAl2024Modelling. <br />

**To simulate the microscopic model:** <br />
- 'MICRO_SimMC_1D.m' : file to run Monte Carlo simulations of the microscopic model in 1D, with Kappa = Dirac delta.  <br />
- 'MICRO_SimMC_1D_tot.m' : file to run Monte Carlo simulations of the microscopic model in 1D, with Kappa = Dirac delta, including interactions in which both phenotypic switching and directional changes occur (results shown in Section 4.2.1). This is an effect of order Dt^2, and in our model set up it yields the same results as in 'MICRO_SimMC_1D.m' (this was thoroughly checked), because phenotypic changes only depend on phenotype and not on velocity, thus simulating the phenotypic switch before the change in velocity accounts for both particles only switching phenotype and those switching both phenotype and velocity. <br />
- 'MICRO_SimMC_2D.m' : file to run Monte Carlo simulations of the microscopic model in 2D (results shown in Section 4.2.1), with Kappa = Dirac delta or Kappa = Von Mises. To run simulations for Kappa given as a Von Mises distributions, first download the function 'vmrand'* by Dylan Muir (2024) and add it to the folder. <br />

* Dylan Muir (2024). vmrand(fMu, fKappa, varargin) (https://www.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin), MATLAB Central File Exchange. Retrieved October 22, 2024.


### Simulating the macroscopic model

The code is set up to simulate the PDE (52) in LorenziEtAl2024Modelling, under the numerical setup detailed in Section 4.1 of LorenziEtAl2024Modelling. 

**To simulate the macroscopic model:** <br />
- 'MACRO_SimPDE_1D.m' : file to run for simulations of the macroscopic model in 1D (results shown in Section 4.2.1) <br />
- 'MACRO_SimPDE_2D.m' : file to run for simulations of the macroscopic model in 2D (results shown in Section 4.2.2) <br />
- 'Nonlocal_advection_1D.m' : function called from within the 'MACRO_Sim_1D.m' file to calculate the advection velocity U_T and the variance-covariance matrix D_T appearing in PDE (52) in 1D <br />
- 'Nonlocal_advection_2D.m' : function called from within the 'MACRO_Sim_2D.m' file to calculate the advection velocity U_T and the variance-covariance matrix D_T appearing in PDE (52) in 2D <br />
- 'MUSCL_GP.m' : function called from within 'MACRO_Sim_1D.m' and 'MACRO_Sim_2D.m' to compute the numerical approximation of the flux of PDE (52) using the MUSCL scheme. Go to this function to change the flux-limiter to use for the simulation. <br />


### Reproducing plots from LorenziEtAl2024Modelling

The data produced by simulations that is used to obtain the plots in LorenziEtAl2024Modelling was stored in a folder 'Saved data'. Due to data file sizes, these are not uploaded in the repository. We detail how to reproduce the figures from scratch and note that the files below to create plots call data files stored in a folder 'Saved data', which may need to be adapted if produced data is stored elsewhere (or with different file names). <br />

**To reproduce the figures in Section 4.2 of LorenziEtAl2024Modelling:** <br />
- 'Plot_comparison_1D_eps.m' : file to reproduce figure 2. One must first run  'MICRO_SimMC_1D_tot.m' (or  'MICRO_SimMC_1D.m') and 'MACRO_SimPDE_1D.m' for the different values of 'epsilon' that should be plot and save the data.
- 'Plot_comparison_2D.m' : file to reproduce figure 3. One must first run  'MICRO_SimMC_2D.m' and 'MACRO_SimPDE_2D.m' for the chosen value of 'epsilon' that should be plot and save the data.
- 'Plot_stripes.m' : function to reproduce figures 4 and 5. This function is called from within 'MACRO_Sim_2D.m' to plot results, so just run that file with the chosen set up. <br />
