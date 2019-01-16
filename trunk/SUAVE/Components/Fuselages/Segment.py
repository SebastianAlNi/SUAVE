## @ingroup Components-Fuselages
# Segment.py
# 
# Created:  Sep 2016, E. Botero (for wings)
# Modified: Jul 2017, M. Clarke
#           Aug 2018, T. St Francis (for fuselages)
#           Jan 2019, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data
from SUAVE.Components import Component, Lofted_Body, Mass_Properties

# ------------------------------------------------------------ 
#  Fuselage Segments
# ------------------------------------------------------------

## @ingroup Components-Fuselages
class Segment(Lofted_Body.Segment):
    def __defaults__(self):
        """This sets the defaults for fuselage segments in SUAVE.
        Assumptions:
        None
        
        Source:
        N/A
        
        Inputs:
        None
        
        Outputs:
        None
        
        Properties Used:
        N/A
        """         
        
        self.tag                = 'segment'
        self.percent_x_location = 0.0      # Percent location along fuselage length.
        self.percent_z_location = 0.0      # Vertical translation of segment. Percent of length.
        self.height             = 0.0
        self.width              = 0.0
        self.length             = 0.0    
        self.effective_diameter = 0.0