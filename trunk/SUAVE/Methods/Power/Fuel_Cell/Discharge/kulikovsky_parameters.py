## @ingroup Methods-Power-Fuel_Cell-Discharge
# kulikovsky_parameters.py
#
# Created : Dec 2019, S. Nicolay
#Adjusted from Matlab model by T. Kadyk (https://www.frontiersin.org/articles/10.3389/fenrg.2019.00035/full) last checked: Dec 2019

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units

# ----------------------------------------------------------------------
#  Setup Larminie
# ----------------------------------------------------------------------

## @ingroup Methods-Power-Fuel_Cell-Discharge
#default values representative of a hydrogen fuel cell
def kulikovsky_parameters(fuel_cell):                     
    """ sets up additional values of fuel cell to run method from Kulikovsky
    
    Inputs:
        fuel cell
     
    Outputs:
        fuel_cell.
          F         [As/mol]
		  l_t       [cm]
		  l_b       [cm]
		  l_m       [cm]
		  b         [V]
		  i_st      [A/cm^3]
		  sigma_t   [S/cm^2]
		  D         [cm^2/s]
		  D_b       [cm^2/s]
		  V_oc      [V]
		  R_o       [Ohm cm^2]
		  c_h       [mol/cm^3]
		  c_ref     [mol/cm^3]
		  
    """   
   
    # natural constants
    fuel_cell.F = 96485.3321  # [As/mol] Faraday constant Ref: NIST CODATA recommendet values
    
    ## geometry
    
    fuel_cell.l_t = 0.001     # [cm] CL thickness    
    fuel_cell.l_b = 0.025     # [cm] GDL thickness
    
    # not used in the model !?
    # fuel_cell.l_m = 0.0025    # [cm] membrane thickness
    
    ## material 
    
    # catalyst
    fuel_cell.b = 0.03        # [V] Tafel slope
    fuel_cell.i_st = 0.817e-3 # [A/cm^3] volumetric exchange currend density
    fuel_cell.sigma_t = 0.03  # [S/cm^2] CCL proton conductivity
    fuel_cell.D = 1.36e-4     # [cm^2/s] oxygen diffusion coefficient in the CL
    
    # GDL
    fuel_cell.D_b = 0.0259    # [cm^2/s] oxygen diffusion coefficient in the GDL
    
    ## cell
    fuel_cell.V_oc = 1.145    # [V] open circuit voltage
    fuel_cell.R_o = 0.126     # [Ohm cm^2] cell ohmic resistivity
    
    ## operation
    fuel_cell.c_h = 7.36e-6# [mol/cm^3] oxygen concentration
    
    ## others
    # the big question is: what's the value for c_ref???
    # own: O2 concentration under standard conditions of air: 101326 Pa, 298K, 21# O2
    # -> this seems to reproduce Fig. 5, curve 1 in the paper
    fuel_cell.c_ref = 8.58335e-6# [mol/cm^3] reference oxygen concentration
    
    # pure O2 under standard conditions:
    # fuel_cell.c_ref = 4.0874e-05# [mol/cm^3] reference oxygen concentration
    
    ## Data from Klingele et al.
    # fuel_cell.V_oc = 0.978489175495164
    # fuel_cell.i_st = 0.042198681176838504
    # fuel_cell.sigma_t = 0.07285302593659942
    # fuel_cell.D = 0.00019398756046993778
    # fuel_cell.D_b = 0.03848596206683334
    # fuel_cell.R_o = 0
    # fuel_cell.c_h = fuel_cell.c_ref
    
    ## Data from Breitwieser et al.
    # fuel_cell.i_st = 0.01221465483410085
    # fuel_cell.sigma_t = 0.005014409221902008
    # fuel_cell.D = 0.0005316516931582585
    # fuel_cell.D_b = 0.03854105682770301
    # fuel_cell.R_o = 0
   
    return