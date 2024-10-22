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
- 'MICRO_SimMC_1D.m' : file to run Monte Carlo simulations of the microscopic model in 1D  <br />
- 'MICRO_SimMC_1D_tot.m' : file to run Monte Carlo simulations of the microscopic model in 1D, including interactions in which both phenotypic switching and directional changes occur (results shown in Section 4.2.1). This is an effect of order Dt^2, and in our model set up it yields the same results as in 'MICRO_SimMC_1D.m' (this was thoroughly checked), because phenotypic changes only depend on phenotype and not on velocity, thus simulating the phenotypic switch before the change in velocity accounts for both particles only switching phenotype and those switching both phenotype and velocity. <br />
- 'MICRO_SimMC_2D.m' : file to run Monte Carlo simulations of the microscopic model in 2D (results shown in Section 4.2.1)  <br />


### Simulating the macroscopic model

The code is set up to simulate the PDE (52) in LorenziEtAl2024Modelling, under the numerical setup detailed in Section 4.1 of LorenziEtAl2024Modelling. 

**To simulate the macroscopic model:** <br />
- 'MACRO_SimPDE_1D.m' : file to run for simulations of the macroscopic model in 1D (results shown in Section 4.2.1) <br />
- 'MACRO_SimPDE_2D.m' : file to run for simulations of the macroscopic model in 2D (results shown in Section 4.2.2) <br />
- 'Nonlocal_advection_1D.m' : function called from within the 'MACRO_Sim_1D.m' file to calculate the advection velocity U_T and the variance-covariance matrix D_T appearing in PDE (52) in 1D <br />
- 'Nonlocal_advection_2D.m' : function called from within the 'MACRO_Sim_2D.m' file to calculate the advection velocity U_T and the variance-covariance matrix D_T appearing in PDE (52) in 2D <br />
- 'MUSCL_GP.m' : function called from within 'MACRO_Sim_1D.m' and 'MACRO_Sim_2D.m' to compute the numerical approximation of the flux of PDE (52) using the MUSCL scheme. Go to this function to change the flux-limiter to use for the simulation. <br />


### Reproducing plots from LorenziEtAl2024Modelling

The data produced by simulations that is used to obtain the plots in LorenziEtAl2024Modelling is provided in the repository in the folder 'Saved data'. We first detail how to reproduce these figures from scratch and then clarify which simulations the different data files are storing results of. <br />

**To reproduce the figures in Section 4.2 of LorenziEtAl2024Modelling:** <br />
- 'Plot_comparison_1D_eps.m' : file to reproduce figure 2. One must first run  'MICRO_SimMC_1D_tot.m' (or  'MICRO_SimMC_1D.m') and 'MACRO_SimPDE_1D.m' for the different values of 'epsilon' that should be plot and save the data.
- 'Plot_comparison_2D.m' : file to reproduce figure 3. One must first run  'MICRO_SimMC_2D.m' and 'MACRO_SimPDE_2D.m' for the chosen value of 'epsilon' that should be plot and save the data.
- 'Plot_stripes.m' : function to reproduce figures 4 and 5. This function is called from within 'MACRO_Sim_2D.m' to plot results, so just run that file with the chosen set up. <br />

**Stored_Data files (.mat):** <br />
- 'datasMC_1D_tot_eps0.0001' : results of  'MICRO_SimMC_1D_tot.m' for epsilon=10^-4 and Kappa='DD' (Figure 2, top row, red) <br />
- 'datasMC_1D_tot_eps0.001' : results of  'MICRO_SimMC_1D_tot.m' for epsilon=10^-3 and Kappa='DD' (Figure 2, middle row, red) <br />
- 'datasMC_1D_tot_eps0.01' : results of  'MICRO_SimMC_1D_tot.m' for epsilon=10^-2 and Kappa='DD' (Figure 2, bottom row, red) <br />
- 'datasMC_2D_eps1e-3_N1e6' : results of  'MICRO_SimMC_2D.m' for epsilon=10^-2 and Kappa='DD' (Figure 3, red) <br />
- 'Saved_11092024_Results_1D_eps0.0001' : results of 'MACRO_SimPDE_1D.m' for epsilon=10^-4 and Kappa='DD' (Figure 2, top row, blue) <br />
- 'Saved_11092024_Results_1D_eps0.001' : results of 'MACRO_SimPDE_1D.m' for epsilon=10^-3 and Kappa='DD' (Figure 2, middle row, blue) <br />
- 'Saved_11092024_Results_1D_eps0.01' : results of 'MACRO_SimPDE_1D.m' for epsilon=10^-2 and Kappa='DD' (Figure 2, bottom row, blue) <br />
- 'Saved_11092024_Setup_1D_eps0.0001' : saved simulation set up yielding  'Saved_11092024_Results_1D_eps0.0001' <br />
- 'Saved_11092024_Setup_1D_eps0.001' : saved simulation set up yielding  'Saved_11092024_Results_1D_eps0.001' <br />
- 'Saved_110924_Results_DD_eps0.001' : results of 'MACRO_SimPDE_2D.m' for epsilon=10^-3 and Kappa='DD' (Figure 4 + Figure 3, blue) <br />
- 'Saved_110924_Results_VM_eps0.001' : results of 'MACRO_SimPDE_2D.m' for epsilon=10^-3 and Kappa='VM' (Figure 5) <br />
- 'Saved_110924_Setup_DD_eps0.001' : saved simulation set up yielding  'Saved_110924_Results_DD_eps0.001' <br />
- 'Saved_110924_Setup_VM_eps0.001' : saved simulation set up yielding  'Saved_110924_Results_VM_eps0.001' <br />