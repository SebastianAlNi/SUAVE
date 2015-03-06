# Tim Momose, January 2015
# Modified February 6, 2015

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pylab as plt
import numpy as np
from copy import deepcopy

# SUAVE Imports
from SUAVE.Core        import Units
from full_setup_737800 import vehicle_setup
from SUAVE.Analyses.Missions.Segments.Conditions import Aerodynamics

# SUAVE-AVL Imports
from SUAVE.Analyses.Aerodynamics.AVL_Callable  import AVL_Callable


def main():

    # -------------------------------------------------------------
    #  Test Script
    # -------------------------------------------------------------

    # Time the process
    import time
    t0 = time.time()
    print "Start: " + time.ctime()

    # Set up test defaults
    vehicle        = vehicle_setup()
    avl            = AVL_Callable()
    avl.keep_files = True
    avl.initialize(vehicle)

        # set up conditions    
    run_conditions = Aerodynamics()
    run_conditions.weights.total_mass[0,0]     = vehicle.mass_properties.max_takeoff
    run_conditions.freestream.mach_number[0,0] = 0.2
    run_conditions.freestream.velocity[0,0]    = 150 * Units.knots
    run_conditions.freestream.density[0,0]     = 1.225
    run_conditions.freestream.gravity[0,0]     = 9.81

    # Set up run cases
    alphas    = np.array([-10,-5,-2,0,2,5,10,20])
    run_conditions.expand_rows(alphas.shape[0])
    run_conditions.aerodynamics.angle_of_attack[:,0] = alphas

    results = avl(run_conditions)

    # Results
    plt.figure('Induced Drag Polar')
    axes = plt.gca()
    CL  = results.aerodynamics.lift_coefficient
    CDi = results.aerodynamics.drag_breakdown.induced.total
    CM  = results.aerodynamics.pitch_moment_coefficient
    axes.plot(CDi,CL,'bo-')
    axes.set_xlabel('Total Drag Coefficient')
    axes.set_ylabel('Total Lift Coefficient')
    axes.grid(True)

    plt.figure('Pitching Momoent')
    axes = plt.gca()
    axes.plot(alphas,CM,'bo-')
    axes.set_xlabel('Angle of Attack')
    axes.set_ylabel('Pitching Moment')
    axes.grid(True)

    tf = time.time()
    print "End:   " + time.ctime()
    print "({0:.2f} seconds)".format(tf-t0)

    plt.show()

    return




# ----------------------------------------------------------------------        
#   Call Main
# ----------------------------------------------------------------------    

if __name__ == '__main__':
    main()