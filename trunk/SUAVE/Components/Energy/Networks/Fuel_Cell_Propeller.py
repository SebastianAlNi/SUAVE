## @ingroup Components-Energy-Networks
# Fuel_Cell_Propeller.py
# 
# Created:  Nov 2019, S. Nicolay
# Derived from Battery:Propeller.py

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import numpy as np
from SUAVE.Components.Propulsors.Propulsor import Propulsor

from SUAVE.Core import Data, Units
from mpi4py import MPI

comm = MPI.COMM_WORLD

# ----------------------------------------------------------------------
#  Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Fuel_Cell_Propeller(Propulsor):
    """ This is a simple network with a fuel cell powering a propeller through
        an electric motor
        
        This network adds 2 extra unknowns to the mission. (The first is
        a voltage, to calculate the thevenin voltage drop in the pack.)
        The second is torque matching between motor and propeller.
    
        Assumptions:
        None
        
        Source:
        None
    """  
    def __defaults__(self):
        """ This sets the default values for the network to function.
    
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
        self.tag = 'fuel_cell_network'     
        self.motor             = None
        self.propeller         = None
        self.esc               = None
        self.avionics          = None
        self.payload           = None
        self.fuel_cell         = None
        self.tank              = None
        self.nacelle_diameter  = None
        self.engine_length     = None
        self.number_of_engines = None
        self.voltage           = None
        self.thrust_angle      = 0.0
    
    # manage process with a driver function
    def evaluate_thrust(self,state):
        """ Calculate thrust given the current state of the vehicle
    
            Assumptions:
            Caps the throttle at 110% and linearly interpolates thrust off that
    
            Source:
            N/A
    
            Inputs:
            state [state()]
    
            Outputs:
            results.thrust_force_vector [newtons]
            results.vehicle_mass_rate   [kg/s]
            conditions.propulsion:
                rpm                  [radians/sec]
                current              [amps]
                fuel_cell_draw       [watts]
                fuel_cell_energy     [joules]
                voltage_open_circuit [volts]
                voltage_under_load   [volts]
                motor_torque         [N-M]
                propeller_torque     [N-M]
    
            Properties Used:
            Defaulted values
        """          
    
        # unpack
        conditions = state.conditions
        numerics   = state.numerics
        motor      = self.motor
        propeller  = self.propeller
        esc        = self.esc
        avionics   = self.avionics
        payload    = self.payload
        fuel_cell  = self.fuel_cell
        num_engines= self.number_of_engines
        #tank       = self.tank
		
        # Step 1 fuel_cell power
        esc.inputs.voltagein = state.unknowns.fuel_cell_voltage_under_load
        #esc.inputs.voltagein = fuel_cell.voltageout
        # Step 2
        esc.voltageout(conditions)
        
        # link
        motor.inputs.voltage = esc.outputs.voltageout
        
        # step 3
        motor.omega(conditions)
        
        # link
        propeller.inputs.omega =  motor.outputs.omega
        propeller.thrust_angle = self.thrust_angle
        
        # step 4
        F, Q, P, Cp, outputs , etap = propeller.spin(conditions)

        #rank = comm.Get_rank()
        #print('Rank: ', rank)
        #print('Cp: ', Cp[0])
        
        # Check to see if magic thrust is needed, the ESC caps throttle at 1.1 already
        eta        = conditions.propulsion.throttle[:,0,None]
        P[eta>1.0] = P[eta>1.0]*eta[eta>1.0]
        F[eta>1.0] = F[eta>1.0]*eta[eta>1.0]
        
        # link
        propeller.outputs = outputs
        
        # Run the avionics
        avionics.power()

        # Run the payload
        payload.power()
        
        # Run the motor for current
        motor.current(conditions)
        
        # link
        esc.inputs.currentout =  motor.outputs.current
        
        # Run the esc
        esc.currentin(conditions)

        # Calculate avionics and payload power
        avionics_payload_power = avionics.outputs.power + payload.outputs.power

        # Calculate avionics and payload current
        avionics_payload_current = avionics_payload_power/self.voltage
        
        # link
        fuel_cell.inputs.current  = esc.outputs.currentin*num_engines + avionics_payload_current
        fuel_cell.inputs.power_in = -(esc.inputs.voltagein*esc.outputs.currentin*num_engines + avionics_payload_power)
        mdot = fuel_cell.energy_calc(conditions, numerics)
        
        # Pack the conditions for outputs
        rpm                  = motor.outputs.omega / Units.rpm
        #current              = esc.outputs.currentin
        current              = fuel_cell.inputs.current
        fuel_cell_draw       = fuel_cell.inputs.power_in
        voltage_under_load   = fuel_cell.voltage_under_load
        cell_voltage         = fuel_cell.cell_voltage
        j0                   = fuel_cell.j0
        efficiency           = fuel_cell.efficiency
        omega_motor          = motor.outputs.omega
        omega_prop           = propeller.inputs.omega
		
		# Calculate thrust coefficient
        Ct = F / (conditions.freestream.density * (rpm/60.0)**2 * (propeller.tip_radius*2)**4) # is /60 correct? Or maybe /60 and 2 pi?

        #print(conditions.freestream.density)
        #print(rpm)
        #print(propeller.prop_attributes.tip_radius)
        
        # Create the outputs 
        conditions.propulsion.rpm                  = rpm
        conditions.propulsion.current              = current
        conditions.propulsion.fuel_cell_draw       = fuel_cell_draw
        conditions.propulsion.motor_voltage        = esc.outputs.voltageout
        conditions.propulsion.voltage_under_load   = voltage_under_load
        conditions.propulsion.cell_voltage         = cell_voltage
        conditions.propulsion.j0                   = j0
        conditions.propulsion.fuel_cell_efficiency = efficiency
        conditions.propulsion.motor_torque         = motor.outputs.torque
        conditions.propulsion.omega_motor          = omega_motor
        conditions.propulsion.propeller_torque     = Q
        conditions.propulsion.propeller_thrust_coefficient = Ct
        conditions.propulsion.propeller_power      = P
        conditions.propulsion.motor_power          = current * esc.outputs.voltageout * conditions.propulsion.etam
        conditions.propulsion.esc_current_in       = esc.outputs.currentin
        conditions.propulsion.esc_current_out      = esc.inputs.currentout
        conditions.propulsion.esc_voltage_in       = esc.inputs.voltagein
        conditions.propulsion.esc_voltage_out      = esc.outputs.voltageout
        
        F    = self.number_of_engines * F * [np.cos(self.thrust_angle),0,-np.sin(self.thrust_angle)]
        conditions.propulsion.mdot = mdot

        results = Data()
        #results.thrust_coefficient  = Ct
        results.thrust_force_vector = F
        results.vehicle_mass_rate   = mdot
        
        return results
    
    
    def unpack_unknowns(self,segment):
        """ This is an extra set of unknowns which are unpacked from the mission solver and send to the network.
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            state.unknowns.propeller_power_coefficient [None]
            state.unknowns.fuel_cell_voltage_under_load  [volts]
    
            Outputs:
            state.conditions.propulsion.propeller_power_coefficient [None]
            state.conditions.propulsion.fuel_cell_voltage_under_load  [volts]
    
            Properties Used:
            N/A
        """                  
        
        # Here we are going to unpack the unknowns (Cp, V) provided for this network
        segment.state.conditions.propulsion.propeller_power_coefficient = segment.state.unknowns.propeller_power_coefficient
        segment.state.conditions.propulsion.fuel_cell_voltage_under_load  = segment.state.unknowns.fuel_cell_voltage_under_load
        
        return
    
    def residuals(self,segment):
        """ This packs the residuals to be send to the mission solver.
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            state.conditions.propulsion:
                motor_torque                          [N-m]
                propeller_torque                      [N-m]
                voltage_under_load                    [volts]
            state.unknowns.fuel_cell_voltage_under_load [volts]
            
            Outputs:
            None
    
            Properties Used:
            self.voltage                              [volts]
        """        
        
        # Here we are going to pack the residuals (torque,voltage) from the network
        
        # Unpack
        q_motor   = segment.state.conditions.propulsion.motor_torque
        q_prop    = segment.state.conditions.propulsion.propeller_torque
        p_prop    = segment.state.conditions.propulsion.propeller_power
        p_motor  = segment.state.conditions.propulsion.motor_power
        v_actual  = segment.state.conditions.propulsion.voltage_under_load
        v_predict = segment.state.unknowns.fuel_cell_voltage_under_load
        v_max     = self.voltage
        
        # Return the residuals
        segment.state.residuals.fuel_cell_network[:,0] = q_motor[:,0] - q_prop[:,0]
        segment.state.residuals.fuel_cell_network[:,1] = (v_predict[:,0] - v_actual[:,0])#/v_max
        
        return    
            
    __call__ = evaluate_thrust
