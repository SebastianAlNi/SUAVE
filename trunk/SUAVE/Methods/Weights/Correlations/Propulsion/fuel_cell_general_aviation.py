## @ingroup Methods-Weights-Correlations-Propulsion
# fuel_cell_general_aviation.py
# 
# Created:  Dec 2019, S. Nicolay
# Modified:

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Units


# ----------------------------------------------------------------------
#   Fuel Cell Propulsion
# ----------------------------------------------------------------------
## @ingroup Methods-Weights-Correlations-Propulsion
def fuel_cell_general_aviation(motor, fuel_cell, tank, bop):
    """ 
        Calculate the weight of the entire fuel cell propulsion system        

        Source:
                nones
                
        Inputs:
                motor mass                                    [kilograms]
                fuel cell mass                                [kilograms]
                tank mass                                     [kilograms]
                balance of plant mass                         [kilograms]
        
        Outputs:
                mass - mass of the full propulsion system     [kilograms]

    """     
    
    mass       = motor + fuel_cell + tank + bop #c all values in kg

    return mass