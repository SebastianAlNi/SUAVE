## @ingroup Methods-Power-Fuel_Cell-Discharge
# kulikovsky.py
#
# Created : Dec 2019, S. Nicolay
# Adjusted from Matlab model by T. Kadyk (https://www.frontiersin.org/articles/10.3389/fenrg.2019.00035/full) last checked: Dec 2019
  
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import scipy as sp
from SUAVE.Core import Units
from .find_voltage_kulikovsky import find_voltage_kulikovsky
from .kulikovsky_parameters import kulikovsky_parameters

# ----------------------------------------------------------------------
#  Kulikovsky
# ----------------------------------------------------------------------

## @ingroup Methods-Power-Fuel_Cell-Discharge
def kulikovsky(fuel_cell,conditions,numerics):
    '''
    function that determines the mass flow rate based on a required power input
    
    Assumptions:
    None (calls other functions)
    
    Inputs:
    fuel_cell.
      area_parallel        [cm^2]
      cells_in_series      [-]
      inputs.
        current            [A]
	  propellant.
	    specific_energy    [J/kg]
    
    Outputs:
	fuel_cell.
	  efficiency           [-]
	  voltage_under_load   [V]
    mdot                   [kg/s]
    
    '''
    
    # Parameters
    p = kulikovsky_parameters(fuel_cell)

    j0 = np.divide(fuel_cell.inputs.current, fuel_cell.area_parallel)
    fuel_cell.j0 = j0

    # Polarization curve
    V_cell = find_voltage_kulikovsky(j0, fuel_cell)
    fuel_cell.cell_voltage = V_cell

    # efficiency
    #fuel_cell.efficiency = np.divide(V_cell, 1.253) # production of water vapor
    # fuel_cell.efficiency_th = V_cell/1.481 % production of liquid water
    fuel_cell.efficiency = np.divide(V_cell, 1.48) # production of liquid water

    #P_nd = np.divide(P, max(P))
    #eta_nd = np.divide(eta, max(eta))
	
    # Pack outputs
    fuel_cell.voltage_under_load   = np.multiply(V_cell, fuel_cell.cells_in_series)
    mdot_fc = np.divide(np.multiply(fuel_cell.voltage_under_load, fuel_cell.inputs.current), np.multiply(fuel_cell.propellant.specific_energy, fuel_cell.efficiency))
    
    mdot    = np.divide(mdot_fc, fuel_cell.fuel_utilization)
   
    return mdot