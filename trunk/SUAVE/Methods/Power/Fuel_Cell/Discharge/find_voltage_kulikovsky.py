## @ingroup Methods-Power-Fuel_Cell-Discharge
# find_voltage_kulikovsky.py
#
# Created : Dec 2019, S. Nicolay
#Adjusted from Matlab model by T. Kadyk (https://www.frontiersin.org/articles/10.3389/fenrg.2019.00035/full) last checked: Dec 2019
  
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from SUAVE.Core import Units

# ----------------------------------------------------------------------
#  Find Voltage Kulikovsky
# ----------------------------------------------------------------------

## @ingroup Methods-Power-Fuel_Cell-Discharge
def find_voltage_kulikovsky(j0, fuel_cell):
    '''
    function that determines the fuel cell voltage based on an input
    current density and some semi-empirical values to describe the voltage
    drop off with current
    
    Assumptions:
    voltage curve is a function of current density of the form
    v = Eoc-r*i1-A1*np.log(i1)-m*np.exp(n*i1)
    
    Inputs:
    current_density           [A/m**2]
    fuel_cell.
        r                     [Ohms*m**2]
        A1                    [V]
        m                     [V]
        n                     [m**2/A]
        Eoc                   [V]
   
    Outputs:
        V                     [V]
         
    
    '''
        #POLCURVE  Analytical Fuel Cell Model after Kulikovsky 2014 JES 161(3) F263-70
    #
    #   [V_cell] = polcurve(j_0,p)  returns the cell voltage V_cell for a given vector of cell current densities j_0   
    #
    #   [V_cell, eta_0] = polcurve(j_0,p)  additionally returns the cell overpotential eta_0
    #
    #   [V_cell, eta_0, P] = polcurve(j_0,p)  additionally returns the cell power density
    #
    
    ## Help Equations
    j_sigma = np.sqrt(2*fuel_cell.i_st*fuel_cell.sigma_t*fuel_cell.b)
    j_st = fuel_cell.sigma_t*fuel_cell.b/fuel_cell.l_t
    beta = np.sqrt(2*j0/j_st)/(1+np.sqrt(1.12*j0/j_st)*np.exp(np.sqrt(2*j0/j_st))) + np.pi*j0/j_st/(2+j0/j_st)
    j_lim = (4*fuel_cell.F*fuel_cell.D_b*fuel_cell.c_h)/fuel_cell.l_b
    
    
    ## Overpotential
    
    eta0 = fuel_cell.b*np.arcsinh(np.power((j0/j_sigma), 2.)/(2*(fuel_cell.c_h/fuel_cell.c_ref)*(1-np.exp(-j0/(2*j_st))))) + fuel_cell.sigma_t*(np.power(fuel_cell.b, 2))/(4.*fuel_cell.F*fuel_cell.D*fuel_cell.c_h)* (j0/j_st - np.log(1+np.power(j0, 2)/(np.power(j_st, 2.)*np.power(beta, 2))))/(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref))) - fuel_cell.b*np.log(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref))) # V
	
    '''if nargout>2:
        eta_actpt = fuel_cell.b*np.arcsinh((j0/j_sigma)^2./(2*(fuel_cell.c_h/fuel_cell.c_ref)*(1-np.exp(-j0/(2*j_st)))))
        eta_tCCL = fuel_cell.sigma_t*fuel_cell.b^2/(4*fuel_cell.F*fuel_cell.D*fuel_cell.c_h)* (j0/j_st - np.log(1+j0^2./(j_st^2*beta^2)))/(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref)))
        eta_tGDL = - fuel_cell.b*np.log(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref)))
        iR = fuel_cell.R_o*j0
        varargout[3] = eta_actpt
        varargout[4] = eta_tCCL
        varargout[5] = eta_tGDL
        varargout[6] = iR
        varargout[7] = beta'''
    eta_actpt = fuel_cell.b*np.arcsinh(np.power((j0/j_sigma), 2.)/(2*(fuel_cell.c_h/fuel_cell.c_ref)*(1-np.exp(-j0/(2*j_st)))))
    eta_tCCL = fuel_cell.sigma_t*np.power(fuel_cell.b, 2)/(4*fuel_cell.F*fuel_cell.D*fuel_cell.c_h)* (j0/j_st - np.log(1+np.power(j0 ,2.)/(np.power(j_st, 2)*np.power(beta, 2))))/(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref)))
    eta_tGDL = - fuel_cell.b*np.log(1- j0/(j_lim*(fuel_cell.c_h/fuel_cell.c_ref)))
    iR = fuel_cell.R_o*j0
    fuel_cell.eta_actpt = eta_actpt
    fuel_cell.eta_tCCL = eta_tCCL
    fuel_cell.eta_tGDL = eta_tGDL
    fuel_cell.iR = iR
    fuel_cell.beta = beta
    
    ## Cell Voltage (V)    
    V_cell = fuel_cell.V_oc - eta0 - fuel_cell.R_o*j0  # [V]
    
    
    ## Power Density (W/cm^2)    
    P = V_cell * j0
    
    
    ## Outputs
    fuel_cell.eta0 = eta0
    fuel_cell.P = P
	
    return V_cell