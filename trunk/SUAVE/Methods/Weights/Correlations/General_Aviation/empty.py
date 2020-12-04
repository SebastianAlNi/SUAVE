## @ingroup Methods-Weights-Correlations-General_Aviation
# empty.py
# 
# Created:  Feb 2018, M. Vegh
# Modified: Dec 2019, S. Nicolay

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Units, Data
from .fuselage import fuselage
from .landing_gear import landing_gear
from .payload import payload
from .systems import systems
from .tail_horizontal import tail_horizontal
from .tail_vertical import tail_vertical
from .wing_main import wing_main
from SUAVE.Methods.Weights.Correlations import Propulsion as Propulsion
import warnings
#from SUAVE.Methods.Center_of_Gravity import compute_aircraft_center_of_gravity

# ----------------------------------------------------------------------
#  Empty
# ----------------------------------------------------------------------
## @ingroup Methods-Weights-Correlations-General_Aviation
def empty(vehicle):
    """ output = SUAVE.Methods.Weights.Correlations.Tube_Wing.empty(engine,wing,aircraft,fuselage,horizontal,vertical)
        Computes the empty weight breakdown of a General Aviation type aircraft  
        
        Inputs:
            engine - a data dictionary with the fields:                    
                thrust_sls - sea level static thrust of a single engine [Newtons]

            vehicle - a data dictionary with the fields:                    
                reference_area                                                            [meters**2]
                envelope - a data dictionary with the fields:
                    ultimate_load - ultimate load of the aircraft                         [dimensionless]
                    limit_load    - limit load factor at zero fuel weight of the aircraft [dimensionless]
                
                mass_properties - a data dictionary with the fields:
                    max_takeoff   - max takeoff weight of the vehicle           [kilograms]
                    max_zero_fuel - maximum zero fuel weight of the aircraft    [kilograms]
                    cargo         - cargo weight                                [kilograms]
                
                passengers - number of passengers on the aircraft               [dimensionless]
                        
                design_dynamic_pressure - dynamic pressure at cruise conditions [Pascal]
                design_mach_number      - mach number at cruise conditions      [dimensionless]
                
                propulsors - a data dictionary with the fields: 
                    keys           - identifier for the type of propulsor; different types have different fields
                        turbofan
                            thrust_sls - sealevel standard thrust                               [Newtons]             
                        internal_combustion
                            rated_power - maximum rated power of the internal combustion engine [Watts]
                        
                    number_engines - integer indicating the number of engines on the aircraft

                wt_cargo - weight of the bulk cargo being carried on the aircraft [kilograms]
                num_seats - number of seats installed on the aircraft [dimensionless]
                ctrl - specifies if the control system is "fully powered", "partially powered", or not powered [dimensionless]
                ac - determines type of instruments, electronics, and operating items based on types: 
                    "short-range", "medium-range", "long-range", "business", "cargo", "commuter", "sst" [dimensionless]
                w2h - tail length (distance from the airplane c.g. to the horizontal tail aerodynamic center) [meters]
                
                fuel - a data dictionary with the fields: 
                    mass_properties  - a data dictionary with the fields:
                        mass -mass of fuel [kilograms]
                    density          - gravimetric density of fuel                             [kilograms/meter**3]    
                    number_of_tanks  - number of external fuel tanks                           [dimensionless]
                    internal_volume  - internal fuel volume contained in the wing              [meters**3]
                wings - a data dictionary with the fields:    
                    wing - a data dictionary with the fields:
                        span                      - span of the wing                           [meters]
                        taper                     - taper ratio of the wing                    [dimensionless]
                        thickness_to_chord        - thickness-to-chord ratio of the wing       [dimensionless]
                        chords - a data dictionary with the fields:
                            mean_aerodynamic - mean aerodynamic chord of the wing              [meters]
                            root             - root chord of the wing                          [meters]
                            
                            
                        sweeps - a data dictionary with the fields:
                            quarter_chord - quarter chord sweep angle of the wing              [radians]
                        mac                       - mean aerodynamic chord of the wing         [meters]
                        r_c                       - wing root chord                            [meters]
                        origin  - location of the leading edge of the wing relative to the front of the fuselage                                      [meters,meters,meters]
                        aerodynamic_center - location of the aerodynamic center of the horizontal_stabilizer relative to the leading edge of the wing [meters,meters,meters]
        
                    
                    
                    
                    horizontal_stabilizer - a data dictionary with the fields:
                        areas -  a data dictionary with the fields:
                            reference - reference area of the horizontal stabilizer                                    [meters**2]
                            exposed  - exposed area for the horizontal tail                                            [meters**2]
                        taper   - taper ratio of the horizontal stabilizer                                             [dimensionless]
                        span    - span of the horizontal tail                                                          [meters]
                        sweeps - a data dictionary with the fields:
                            quarter_chord - quarter chord sweep angle of the horizontal stabilizer                     [radians]
                        chords - a data dictionary with the fields:
                            mean_aerodynamic - mean aerodynamic chord of the horizontal stabilizer                     [meters]         
                            root             - root chord of the horizontal stabilizer             
                        thickness_to_chord - thickness-to-chord ratio of the horizontal tail                           [dimensionless]
                        mac     - mean aerodynamic chord of the horizontal tail                                        [meters]
                        origin  - location of the leading of the horizontal tail relative to the front of the fuselage                                                 [meters,meters,meters]
                        aerodynamic_center - location of the aerodynamic center of the horizontal_stabilizer relative to the leading edge of the horizontal stabilizer [meters,meters,meters]
        
                    vertical - a data dictionary with the fields:
                        areas -  a data dictionary with the fields:
                            reference - reference area of the vertical stabilizer         [meters**2]
                        span    - span of the vertical tail                               [meters]
                        taper   - taper ratio of the horizontal stabilizer                [dimensionless]
                        t_c     - thickness-to-chord ratio of the vertical tail           [dimensionless]
                        sweeps   - a data dictionary with the fields:
                            quarter_chord - quarter chord sweep angle of the vertical stabilizer [radians]
                        t_tail - flag to determine if aircraft has a t-tail, "yes"               [dimensionless]


                
                fuselages - a data dictionary with the fields:  
                    fuselage - a data dictionary with the fields:
                        areas             - a data dictionary with the fields:
                            wetted - wetted area of the fuselage [meters**2]
                        differential_pressure  - Maximum fuselage pressure differential   [Pascal]
                        width             - width of the fuselage                         [meters]
                        heights - a data dictionary with the fields:
                            maximum - height of the fuselage                              [meters]
                        lengths-  a data dictionary with the fields:
                            structure - structural length of the fuselage                 [meters]                     
                        mass_properties - a data dictionary with the fields:
                            volume - total volume of the fuselage                         [meters**3]
                            internal_volume - internal volume of the fuselage             [meters**3]
                        number_coach_sets - number of seats on the aircraft               [dimensionless]    
                landing_gear - a data dictionary with the fields:
                    main - a data dictionary with the fields:
                        strut_length - strut length of the main gear                      [meters]
                    nose - a data dictionary with the fields:
                        strut_length - strut length of the nose gear                      [meters]
                avionics - a data dictionary, used to determine if avionics weight is calculated, don't include if vehicle has none
                air_conditioner - a data dictionary, used to determine if air conditioner weight is calculated, don't include if vehicle has none
        
        
        Outputs:
            output - a data dictionary with fields:
                wing - wing weight                            [kilograms]
                fuselage - fuselage weight                    [kilograms]
                propulsion - propulsion                       [kilograms]
                landing_gear_main - main gear weight          [kilograms]
                landing_gear_nose - nose gear weight          [kilograms]
                horizonal_tail - horizontal stabilizer weight [kilograms]
                vertical_tail - vertical stabilizer weight    [kilograms]
                systems - total systems weight                [kilograms]
                systems_breakdown - a data dictionary with fields:
                    control_systems - control systems weight  [kilograms]
                    hydraulics - hydraulics weight            [kilograms]
                    avionics - avionics weight                [kilograms]
                    electric - electrical systems weight      [kilograms]
                    air_conditioner - air conditioner weight  [kilograms]
                    furnish - furnishing weight               [kilograms]
                    fuel_system - fuel system weight          [ kilograms]
           Wing, empannage, fuselage, propulsion and individual systems masses updated with their calculated values
       Assumptions:
            calculated aircraft weight from correlations created per component of historical aircraft
        
    """     

    # Unpack inputs
    S_gross_w   = vehicle.reference_area
    fuel        = vehicle.fuel
    Nult        = vehicle.envelope.ultimate_load
    Nlim        = vehicle.envelope.limit_load
    TOW         = vehicle.mass_properties.max_takeoff
    wt_zf       = vehicle.mass_properties.max_zero_fuel
    num_pax     = vehicle.passengers
    wt_cargo    = vehicle.mass_properties.cargo

    q_c         = vehicle.design_dynamic_pressure
    mach_number = vehicle.design_mach_number

    propulsor_name = list(vehicle.propulsors.keys())[0] #obtain the key for the propulsor for assignment purposes
    propulsors     = vehicle.propulsors[propulsor_name]
    num_eng        = propulsors.number_of_engines

    if propulsor_name == 'turbofan':
        thrust_sls                       = propulsors.sealevel_static_thrust
        wt_engine_jet                    = Propulsion.engine_jet(thrust_sls)
        wt_propulsion                    = Propulsion.integrated_propulsion(wt_engine_jet,num_eng)
        propulsors.mass_properties.mass  = wt_propulsion 

    elif propulsor_name == 'internal_combustion':
        rated_power                      = propulsors.rated_power/num_eng
        wt_engine_piston                 = Propulsion.engine_piston(rated_power)
        wt_propulsion                    = Propulsion.integrated_propulsion_general_aviation(wt_engine_piston,num_eng)
        propulsors.mass_properties.mass  = wt_propulsion 
		
    #added by S.Nicolay
    elif propulsor_name == 'fuel_cell_network':
        wt_motor                         = propulsors.motor.mass_properties.mass * num_eng
        wt_fuel_cell                     = propulsors.fuel_cell.mass_properties.mass
        wt_tank                          = propulsors.tank.mass_properties.mass
        wt_bop                           = propulsors.bop.mass_properties.mass
        wt_propulsion                    = Propulsion.fuel_cell_general_aviation(wt_motor, wt_fuel_cell, wt_tank, wt_bop)
        propulsors.mass_properties.mass  = wt_propulsion

    else: #propulsor used is not an IC Engine or Turbofan; assume mass_properties defined outside model
        wt_propulsion                    = propulsors.mass_properties.mass
        if wt_propulsion==0:
            warnings.warn("Propulsion mass= 0 ;e there is no Engine Weight being added to the Configuration", stacklevel=1)    
    #find fuel volume
    if 'fuel' not in vehicle: 
        warnings.warn("fuel mass= 0 ; fuel system volume is calculated incorrectly ", stacklevel=1)   
        N_tank     = 0 
        V_fuel     = 0.
        V_fuel_int = 0.
      
    else:
        m_fuel                      = fuel.mass_properties.mass
        landing_weight              = TOW-m_fuel  #just assume this for now
        N_tank                      = fuel.number_of_tanks
        V_fuel_int                  = fuel.internal_volume #fuel in internal (as opposed to external tanks)
        V_fuel                      = m_fuel/fuel.density  #total fuel
        fuel.mass_properties.volume = V_fuel 

    #main wing
    if 'main_wing' not in vehicle.wings:
        wt_wing = 0.0
        wing_c_r = 0.0
        warnings.warn("There is no Wing Weight being added to the Configuration", stacklevel=1)

    else:
        b          = vehicle.wings['main_wing'].spans.projected
        AR_w       = (b**2.)/S_gross_w
        taper_w    = vehicle.wings['main_wing'].taper
        t_c_w      = vehicle.wings['main_wing'].thickness_to_chord
        sweep_w    = vehicle.wings['main_wing'].sweeps.quarter_chord
        mac_w      = vehicle.wings['main_wing'].chords.mean_aerodynamic
        wing_c_r   = vehicle.wings['main_wing'].chords.root
        #now run weight script for the wing
        #wt_wing                                         = wing_main(S_gross_w, m_fuel, AR_w, sweep_w, q_c, taper_w, t_c_w,Nult,TOW)
        # Fuel cell is in the main wing, so use fuel cell mass instead of fuel mass for wing mass estimation
        wt_wing                                         = wing_main(S_gross_w, wt_fuel_cell, AR_w, sweep_w, q_c, taper_w, t_c_w,Nult,TOW)
        vehicle.wings['main_wing'].mass_properties.mass = wt_wing
        
        '''for i in range(10, 41):
            wt_wing     = wing_main(S_gross_w, wt_fuel_cell, i, sweep_w, q_c, taper_w, t_c_w,Nult,TOW)
            print(i, wt_wing.item(0))'''
   
    if 'horizontal_stabilizer' not in vehicle.wings:
        wt_tail_horizontal = 0.0
        S_h = 0.0
        warnings.warn("There is no Horizontal Tail Weight being added to the Configuration", stacklevel=1)

    else:    
        S_h                = vehicle.wings['horizontal_stabilizer'].areas.reference
        b_h                = vehicle.wings['horizontal_stabilizer'].spans.projected
        AR_h               = (b_h**2.)/S_h
        taper_h            = vehicle.wings['horizontal_stabilizer'].spans.projected
        sweep_h            = vehicle.wings['horizontal_stabilizer'].sweeps.quarter_chord
        mac_h              = vehicle.wings['horizontal_stabilizer'].chords.mean_aerodynamic
        t_c_h              = vehicle.wings['horizontal_stabilizer'].thickness_to_chord
        l_w2h              = vehicle.wings['horizontal_stabilizer'].origin[0][0] + vehicle.wings['horizontal_stabilizer'].aerodynamic_center[0] - vehicle.wings['main_wing'].origin[0][0] - vehicle.wings['main_wing'].aerodynamic_center[0] #used for fuselage weight
        wt_tail_horizontal = tail_horizontal(S_h, AR_h, sweep_h, q_c, taper_h, t_c_h,Nult,TOW)                
        
        vehicle.wings['horizontal_stabilizer'].mass_properties.mass = wt_tail_horizontal        
    
    #vertical stabilizer
    if 'vertical_stabilizer' not in vehicle.wings:   
        output_3 = Data()
        output_3.wt_tail_vertical = 0.0

        S_v = 0.0
        warnings.warn("There is no Vertical Tail Weight being added to the Configuration", stacklevel=1)    

    else:     
        S_v        = vehicle.wings['vertical_stabilizer'].areas.reference
        b_v        = vehicle.wings['vertical_stabilizer'].spans.projected
        AR_v       = (b_v**2.)/S_v
        taper_v    = vehicle.wings['vertical_stabilizer'].taper
        t_c_v      = vehicle.wings['vertical_stabilizer'].thickness_to_chord
        sweep_v    = vehicle.wings['vertical_stabilizer'].sweeps.quarter_chord
        t_tail     = vehicle.wings['vertical_stabilizer'].t_tail  
        output_3   = tail_vertical(S_v, AR_v, sweep_v, q_c, taper_v, t_c_v, Nult,TOW,t_tail)
        
        vehicle.wings['vertical_stabilizer'].mass_properties.mass = output_3.wt_tail_vertical
    
    if 'fuselages.fuselage' not in vehicle:
        S_fus      = vehicle.fuselages['fuselage'].areas.wetted
        diff_p_fus = vehicle.fuselages['fuselage'].differential_pressure
        w_fus      = vehicle.fuselages['fuselage'].width
        h_fus      = vehicle.fuselages['fuselage'].heights.maximum
        l_fus      = vehicle.fuselages['fuselage'].lengths.structure
        V_fuse     = vehicle.fuselages['fuselage'].mass_properties.volume
        #V_int      = vehicle.fuselages['fuselage'].mass_properties.internal_volume 
        num_seats  = vehicle.fuselages['fuselage'].number_coach_seats 
        #calculate fuselage weight
        wt_fuselage = fuselage(S_fus, Nult, TOW, w_fus, h_fus, l_fus, l_w2h, q_c, V_fuse, diff_p_fus)
    else:
        print('got here')
        warnings.warn('There is no Fuselage weight being added to the vehicle', stacklevel=1)

    #landing gear
    if 'landing_gear' not in vehicle:
        warnings.warn('There is no Landing Gear weight being added to the vehicle', stacklevel=1)
        wt_landing_gear = Data()
        wt_landing_gear.main = 0.0
        wt_landing_gear.nose = 0.0

    else: 
        landing_gear_component = vehicle.landing_gear #landing gear previously defined
        strut_length_main      = landing_gear_component.main_strut_length
        strut_length_nose      = landing_gear_component.nose_strut_length
        wt_landing_gear        = landing_gear(landing_weight, Nult, strut_length_main, strut_length_nose)
        
        #landing_gear_component.main.mass_properties.mass = wt_landing_gear.main
        #landing_gear_component.nose.mass_properties.mass = wt_landing_gear.nose

    if 'avionics' not in vehicle:
        warnings.warn('There is no avionics weight being added to the vehicle; many weight correlations are dependant on this', stacklevel=1)
        avionics          = SUAVE.Components.Energy.Peripherals.Avionics()
        W_uav = 0.

    else:
        avionics = vehicle.avionics
        W_uav = avionics.mass_properties.uninstalled

    has_air_conditioner = 0
    if 'air_conditioner' in vehicle:
        has_air_conditioner = 1

    # Calculating Empty Weight of Aircraft
    output_2           = systems(W_uav,V_fuel, V_fuel_int, N_tank, num_eng, l_fus, b, TOW, Nult, num_seats, mach_number, has_air_conditioner)

    # Calculate the equipment empty weight of the aircraft

    wt_empty           = (wt_wing + wt_fuselage + wt_landing_gear.main+wt_landing_gear.nose + wt_propulsion + output_2.wt_systems + \
                          wt_tail_horizontal + output_3.wt_tail_vertical) 
    vehicle.fuselages['fuselage'].mass_properties.mass = wt_fuselage

    # packup outputs
    output                   = payload(TOW, wt_empty, num_pax,wt_cargo)
    output.wing              = wt_wing
    output.fuselage          = wt_fuselage
    output.propulsion        = wt_propulsion
    output.landing_gear_main = wt_landing_gear.main
    output.landing_gear_nose = wt_landing_gear.nose
    output.horizontal_tail   = wt_tail_horizontal
    output.vertical_tail     = output_3.wt_tail_vertical

    output.systems           = output_2.wt_systems       
    output.systems_breakdown = Data()
    output.systems_breakdown.control_systems   = output_2.wt_flt_ctrl        
    output.systems_breakdown.hydraulics        = output_2.wt_hyd_pnu     
    output.systems_breakdown.avionics          = output_2.wt_avionics    
    output.systems_breakdown.electrical        = output_2.wt_elec        
    output.systems_breakdown.air_conditioner   = output_2.wt_ac          
    output.systems_breakdown.furnish           = output_2.wt_furnish    
    output.systems_breakdown.fuel_system       = output_2.wt_fuel_sys
	
	# added by S. Nicolay
    output.propulsion        = wt_propulsion
    output.propulsion_breakdown = Data()
    output.propulsion_breakdown.motors         = wt_motor
    output.propulsion_breakdown.fuel_cell      = wt_fuel_cell
    output.propulsion_breakdown.tank           = wt_tank
    output.propulsion_breakdown.bop            = wt_bop	

    #define weights components

    control_systems   = SUAVE.Components.Physical_Component()
    electrical_systems= SUAVE.Components.Physical_Component()
    passengers        = SUAVE.Components.Physical_Component()
    furnishings       = SUAVE.Components.Physical_Component()
  
    apu               = SUAVE.Components.Physical_Component()
    hydraulics        = SUAVE.Components.Physical_Component()

    #assign output weights to objects
    control_systems.mass_properties.mass                             = output.systems_breakdown.control_systems
    electrical_systems.mass_properties.mass                          = output.systems_breakdown.electrical
    passengers.mass_properties.mass                                  = output.pax + output.bag
    furnishings.mass_properties.mass                                 = output.systems_breakdown.furnish
    avionics.mass_properties.mass                                    = output.systems_breakdown.avionics              
    hydraulics.mass_properties.mass                                  = output.systems_breakdown.hydraulics

    if has_air_conditioner:
        vehicle.air_conditioner.mass_properties.mass                 = output.systems_breakdown.air_conditioner

    #assign components to vehicle
    vehicle.control_systems                     = control_systems
    vehicle.electrical_systems                  = electrical_systems
    vehicle.avionics                            = avionics
    vehicle.furnishings                         = furnishings
    vehicle.passenger_weights                   = passengers 
   
    vehicle.apu                                 = apu
    vehicle.hydraulics                          = hydraulics
    vehicle.landing_gear                        = landing_gear_component
	
	# added
    #compute_aircraft_center_of_gravity(vehicle)

    #note; air conditioner optional, and weight is added to the air_conditioner object directly
    return output
