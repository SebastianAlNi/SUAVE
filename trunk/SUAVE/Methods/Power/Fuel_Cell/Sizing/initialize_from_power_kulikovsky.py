## @ingroup Methods-Power-Fuel_Cell-Sizing
# initialize_from_power_kulikovsky.py
#
# Created : Dec. 2019, S. Nicolay

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import scipy as sp
import numpy as np
from SUAVE.Core import Units
           
# ----------------------------------------------------------------------
#  Initialize from Power
# ----------------------------------------------------------------------

## @ingroup Methods-Power-Fuel_Cell-Sizing
def initialize_from_power_kulikovsky(fuel_cell):
    '''
    assigns the mass of the fuel cell based on the power and specific power
    Assumptions:
    specific power and mass densities include both the stack and the periphery
    
    Inputs:
    fuel_cell.
      max_power          [W]
      max_power_density  [W/cm^2]
      specific_power     [W/kg]
      cells_in_series    [-]
      mass_density       [kg/m^3]
    
    Outputs:
    fuel_cell.
      cell_area      [cm^2]
      area_parallel  [cm^2]
      volume         [m^3]
      mass_properties.
        mass         [kg]
    '''
	
    fuel_cell.cell_area            = fuel_cell.max_power/fuel_cell.max_power_density
    fuel_cell.area_parallel        = fuel_cell.cell_area/fuel_cell.cells_in_series
    fuel_cell.mass_properties.mass = fuel_cell.max_power/fuel_cell.specific_power *Units.kg
    fuel_cell.volume               = fuel_cell.mass_properties.mass/fuel_cell.mass_density
