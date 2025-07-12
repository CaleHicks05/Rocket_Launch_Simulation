import numpy as np
import matplotlib.pyplot as plt

# Parameters
m_stage1_full = 430000.0 # Mass of Stage 1 with full fuel
m_stage1_empty = 25000.0 # Mass of Stage 1 w/o fuel
m_stage2_full = 111000.0 # Mass of Stage 2 with full fuel
m_stage2_empty = 4000.0 # Mass of Stage 2 w/o fuel

stage1_burn_time = 150.0 # seconds
stage2_burn_time = 360.0 # seconds

ve1 = 2800 # exhaust velocity of Stage 1
ve2 = 3400 # exhaust velocity of Stage 2

initial_mass = m_stage1_full + m_stage2_full # Total initial mass of rocket
after_stage1_mass = m_stage1_empty + m_stage2_full # Mass w/ empty stage 1

thrust_1 = ve1 * (initial_mass - m_stage1_empty) / stage1_burn_time # Thrust produced by stage 1
thrust_2 = ve2 * (m_stage2_full - m_stage2_empty) / stage2_burn_time # Thrust produced by stage 2

# Constants
g = 9.81 # gravity [m/s^2]
G = 6.67 * (10 ** (-11)) # Gravitational Constant
#Cd = 0.3 # drag coefficient (will change, this is a guess)
# rho = 1.225 # air density [kg/m^3]
A = 10.2 # cross-sectional area of the rocket [m^2]
ve = 50000 # exhaust velocity [m/s]
M_earth = 5.97 * 10 ** 24 # [kg] - mass of earth
R_earth = 6371000 # meters - radius of earth

# Earth Atmospheric Model (NASA)
def atmospheric_model(h): # defines function of NASA's armospheric model
    if h > 25000: # Upper Stratosphere
        T = -131.21 + 0.00299 * h # Temperature in Celsius
        P = 2.488 * ((T + 273.15)/216.6)**(-11.388) # Pressure in kPa
    elif 11000 <= h <= 25000: # Lower Stratosphere
        T = -56.46 # Temperature in Celsius
        P = 22.65 * np.exp(1.73 - 0.000157 * h) # Pressure in kPa
    else: # Troposphere
        T = 15.04 - 0.00649 * h # Tempeature in Celsius
        P = 101.29 * ((T + 273.15)/288.08)**(5.256) # Pressure in kPa
    return T, P # Return temp. and pressure

def air_density(altitude): # defines function of air density using ideal gas law
    T_C, P_kPa = atmospheric_model(altitude) # Call temp (in celsius) and pressure (in kPa) from armospheric function
    T_K = T_C + 273.15 # convert to kelvin
    P_Pa = P_kPa * 1000 # convert to pascals
    R_specific = 287.05 # specific gas constant for air approx.
    return P_Pa / (R_specific * T_K) # return rho based on gas law

def speed_of_sound(T): # defines function of the speed of sound based on the temperature of medium
    gamma = 1.4 # ratio of specific heats between medium for air
    R_specific = 287.05 # specifc gas constant for air
    return (gamma * R_specific * T) ** 0.5 # return the speed of sound depending on temp.

def drag_coefficient(mach): # defines function of drag coefficient using table found
    if mach < 0.5: # low subsonic
        return 0.2
    elif 0.5 <= mach < 0.8: # moderate subsonic
        return 0.3
    elif 0.8 <= mach < 1.2: # transonic
        return 0.6
    elif 1.2 <= mach < 2: # low supersonic
        return 0.4
    elif 2 <= mach < 5: # supersonic
        return 0.35
    elif 5 <= mach < 10: # high supersonic
        return 0.3
    else: # hypersonic
        return 0.25
    
def F_drag(v, h): # defins function of drag force against rocket based on velocity and altitude
    T_C, P_kPa = atmospheric_model(h) # defines h as altitude call-back and temp and pressure, temp is used for drag but pressure is not needed
    T_K = T_C + 273.15 # convert to kelvin
    v_sound = speed_of_sound(T_K) # speed of sound based on temp.
    mach = abs(v) / v_sound # mach number
    Cd = drag_coefficient(mach) # coeficcient of drag based on the speed of the rocket call-back to function
    rho = air_density(h) # call-back to function using h as altitude to find air density
    return 0.5 * rho * A * Cd * v**2 * np.sign(v) # returns drag forc calculation

# Time setup
dt = 0.1 # small change in time [s]
t_max = 50000.0 # total simulation time [s]
steps = int(t_max / dt) # number of simulations/calculations ran

# Initialize arrays
t = np.zeros(steps) # initialize all time arrays to zeros for number of steps ran
altitude = np.zeros(steps) # initialize all altitude arrays to zeros for number of steps ran
velocity = np.zeros(steps) # initialize all velocity arrays to zeros for number of steps ran
mass = np.full(steps, initial_mass) # initialize all mass arrays to the initial mass for the number of steps ran
F_thrust_arr = np.zeros(steps) # initialize all thrust arrays to zeros for number of steps ran
x_positions = [] # initialize horizontal positioning as an empty list to collect x-coords
y_positions = [] # initialize vertical positioning as an empty list to collect y-coords

# parachute_deployment = False # parachute is not deployed

def F_gravity(m, h): # function defins the force of gravity based on mass and altitude
        return G * M_earth * m / (R_earth + h) ** 2 # returns force towards Earth

# Simulation
stage1_ignited = False # stage 1 has not been ignited
stage2_ignited = False # stage 2 has not been ignited
stage1_seperated = False # stage 1 has not seperated yet
orbit_sim_ran = False # orbit simulation has not been ran
apogee_detected = False # apogee has not been reached/detected

for i in range(1, steps): # for loop starts at 1 so time change doesn't go negative and goes until number of steps
    # use velocity[i-1] to calculate forces and acceleration
    # update velocity[i] based on previous velocity and acceleration

    t[i] = t[i-1] + dt # time at i is equal to time at i-1 plus the change in time (dt)

    if t[i] <= stage1_burn_time: # Stage 1 burn if time is less than or equal to its burn time (150 s)
        if not stage1_ignited: # if stage 1 has not been ignited yet
            stage1_ignited = True # stage 1 has begun ignition
        stage = 1 # indexing easier
        mass[i] = initial_mass - (m_stage1_full - m_stage1_empty) * (t[i] / stage1_burn_time) # mass at i is equal to the initial mass minus mass lost over time frame for initial burn of stage 1 only
        F_thrust = thrust_1 # force of thrust is equal to the thrust produced by stage 1 for this time frame
        F_thrust_arr[i] = F_thrust # force array at i is equal to thrust at i
    elif stage1_burn_time < t[i] <= stage1_burn_time + stage2_burn_time: # Stage 2 burn if time is greater than stage 1 time and less than their total combined time
        if not stage2_ignited: # if stage 2 has not been ignited yet
            stage2_ignited = True # stage 2 has begun ignition
        if not stage1_seperated: # if stage 1 has not begun seperation yet
            print(f"Stage 1 seperated at t = {t[i]} s and altitude = {altitude[i-1]} meters") # prints when stafe 1 seperates
            stage1_seperated = True # stage 1 has begun seperation
            mass[i] -= m_stage1_empty # Drops empty first stage
            velocity[i] += 50 # Recoil velocity kick estimate
        stage = 2 # index easier
        stage2_start_time = t[i] - stage1_burn_time # stage 2 elapsed time is equal to total time minus stage 1 time
        mass[i] = m_stage2_full - (m_stage2_full - m_stage2_empty) * (stage2_start_time / stage2_burn_time) # mass at time i is equal to the initial mass of stage 2 full minus mass lost over time fram for stage 2
        F_thrust = thrust_2 # force of thrust is equal to the thrust produced by stage 2 for this time frame
        F_thrust_arr[i] = F_thrust # force array at i is equal to thrust at i
    else: # No burn/after both are done burning
        stage = 0 # momentum with no thrust
        mass[i] = m_stage2_empty # mass is left to be only of empty stage 2
        F_thrust = 0 # no additional thrust
        F_thrust_arr[i] = F_thrust # force array at i is equal to 0

    # Forces
    F_d = F_drag(velocity[i-1], altitude[i-1]) # force of drag call-back to approx. velocity and altitude and drag at time i using i-1 
    F_g = F_gravity(mass[i], altitude[i]) # force of gravity call-back using mass and altitude at time i
    F_net = F_thrust - F_g - F_d # net total force is thrust minus drag minus gravity
    a = F_net / mass[i] # acceleration of rocket is net force / mass at time i
    velocity[i] = velocity[i-1] + a * dt # velocity at time i is equal to previous velocity plus acceleration by change in time (dt)
    altitude[i] = altitude[i-1] + velocity[i] * dt # altitude at time i is equal to previous altitude plus current velocity by change in time (dt)

    #if i % 200 == 0:
        #print(f"Step {i}: altitude = {altitude[i]:.1f} m, velocity = {velocity[i]:.1f} m/s")

# Apogee Detection
    if velocity[i-1] > 0 and velocity[i] <= 0: # if the previous velocity is postive and current velocity is zero or negative it is concave down meaning there is apogee
        apogee_index = i # index apogee at step i
        apogee_time = t[i] # time at apogee is equal to time at i
        apogee_altitude = altitude[i] # altitude at apogee is equal to altitude at i
        apogee_detected = True # apogee has been detected
        print(f"Apogee at altitude = {apogee_altitude} and time = {apogee_time}") # prints time and altitude where apogee occurs

# print(f"Debug: Apogee detected? {apogee_detected}, Orbit sim already ran? {orbit_sim_ran}") # debug statement to find apogee
# Orbit/Tangental Velocity
if apogee_detected and not orbit_sim_ran: # if apogee has been detected and the orbit sim has not been ran
    orbit_sim_ran = True # then orbit sim has begun
    r_min = np.inf  # start with a very large number for distance
    r_max = 0       # start with a very small number
    position = np.array([0.0 , R_earth + apogee_altitude]) # (horizontal, vertical) vector, ignoring z-comp. for position
    v_orbit = np.sqrt((G * M_earth) / (R_earth + apogee_altitude)) # Combination between F_gravity = F_centripetal
    v_centripetal = np.array([v_orbit, 0.0]) # set orbit velocity array to x- v_orbit and y to 0

    t_total = 12000 # time for sim of orbit
    num_steps = int(t_total / dt) # step number for orbit

    for step in range(num_steps): # for step in range of the number of steps for orbit
        r = np.linalg.norm(position) # finds radius from core of earth using euclidian magnitudes
        a_centripetal = (-G * M_earth / r**3 )* position # centripetal acceleration is found using gravitational force and centripetal force
        v_centripetal += a_centripetal * dt # centripetal velocity found using acceleration times change in time (added continuously)
        position += v_centripetal * dt # position in x-coord. found using change in velocity by change in time (added continously)
        orbit_altitude = r - R_earth # height from surface of Earth 
        speed = np.linalg.norm(v_centripetal) # scalar magnitude of vector

        x_positions.append(position[0]) # grab x_positions found and put them into x array
        y_positions.append(position[1]) # grab y-positions found and put them into y array

        # Energy
        U = -G * M_earth * m_stage2_empty / r # potential energy due to gravity of empty stage 2
        K = 0.5 * m_stage2_empty * (speed ** 2) # kinetic energy due to magnitude (speed) of velocity and empty stage 2 rocket
        E_mechanical = U + K # mechanical energy equal to potential plus kinetic energies

        # Track periapsis and apoapsis
        if r < r_min: # if distance is less than min. radius, then it is continuously ran until smallest r is found
            r_min = r # set r equal to new r_min
        if r > r_max: # if distance is greater than max radius, then it is continuously ran until largest r is found
            r_max = r # set r equal to new r_max

        if step % 500 == 0: # for intervals of step divisible by 500 evenly
            print(f"Total mechanical energy is {E_mechanical} J, Orbiting at altitude = {altitude[step]} m, speed = {speed}, at time = {step * dt} s") # print total energy, height, speed, and time

    # Orbital elements
    periapsis = r_min - R_earth # closest point to surface of earth in orbit
    apoapsis = r_max - R_earth # furtherst point to sruface of earth in orbit
    eccentricity = (r_max - r_min) / (r_max + r_min) # curve determined (0 means it is circular)

    print("\n----- Orbital Parameters -----") # print heading
    print(f"Apoapsis: {apoapsis} m") # print apoapsis in meters
    print(f"Periapsis: {periapsis} m") # print periapsis in meters
    print(f"Eccentricity: {eccentricity}") # print eccentricity

if x_positions and y_positions: # if both x and y positions are found 
    plt.figure(figsize=(8, 8)) # graph is square
    plt.plot(x_positions, y_positions, label='Orbit Path') # path of orbit of rocket
    plt.plot(0, 0, 'yo', label='Earth Center') # yellow center dot (orgion)
    earth = plt.Circle((0, 0), R_earth, color='blue', alpha=0.3, label='Earth Surface') # radius is 0.3 of surface of earth
    plt.gca().add_patch(earth) # fill in area 
    plt.xlabel('X Position (m)') # x-label
    plt.ylabel('Y Position (m)') #y-label
    plt.title('Rocket Orbit Trajectory') # title
    plt.axis('equal') # axes are equal in length (sqaure)
    plt.grid(True) # plot grid
    plt.legend() # plot legend
    plt.show() 
else: # if x and y positions cannot be obtained
    print("Orbit path not simulated — skipping plot.") # print orbit cannot be simulated


# Plot Sim
if 1 == 1 is True:
    plt.figure(figsize=(10,6)) # mass depletion
    plt.plot(t, mass, label='Mass (kg)')
    plt.xlabel('Time (s)')
    plt.ylabel('Mass (kg)')
    plt.title('Rocket mass vs Time')
    plt.grid(True)
    plt.legend()
    plt.show

    plt.figure(figsize=(10,6)) # altitude
    plt.plot(t, altitude, label='Altitude (m)')
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    plt.title('Rocket Altitude vs Time')
    plt.grid(True)
    plt.legend()

    plt.figure(figsize=(10,6)) # velocity
    plt.plot(t, velocity, label='Velocity (m/s)', color='orange')
    plt.axhline(0, color='red', linestyle='--', label='Zero Velocity')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Rocket Velocity vs Time')
    plt.grid(True)
    plt.legend()

    plt.figure(figsize=(10, 6)) # thrust
    plt.plot(t, F_thrust_arr, label='Thrust (N)', color='orange')
    plt.xlabel('Time (s)')
    plt.ylabel('Thrust (N)')
    plt.title('Rocket Thrust vs Time')
    plt.grid(True)
    plt.legend()
    plt.show()

    plt.show()
    plt.show(block=True)


#matplotlib.animation


################# Not needed now

# def air_density(altitude_m):
   # rho_0 = 1.225 # sea level
   # H = 1000000 # Scale height
   # return rho_0 *np.exp(-altitude_m / H)
# def F_drag(v, h):
    #rho = air_density(h)
    #return 0.5 * rho * A * Cd * v**2 * np.sign(v)
# Simulation
# for i in range(1, steps):
    # use velocity[i-1] to calculate forces and acceleration
    # update velocity[i] based on previous velocity and acceleration
   # t[i] = t[i-1] + dt
   # if t[i] <= burn_time:
   #     mass[i] = m0 - (m0 - mf) * (t[i] / burn_time)
   #     F_thrust = thrust
   # else: 
   #     mass[i] = mf
   #     F_thrust = 0

   # F_d = F_drag(velocity[i-1], altitude[i-1])
   # F_net = F_thrust - (mass[i] * g) - F_d
   # a = F_net / mass[i]
   # velocity[i] = velocity[i-1] + a * dt
   # altitude[i] = altitude[i-1] + velocity[i] * dt

   # if velocity[i-1] > 0 and velocity[i] <= 0:
   #     apogee_index = 1
   #     apogee_time = t[i]
   #     apogee_altitude = altitude[i]
   #     parachute_deployment = True
   #     print(f"Parachute deployed: {parachute_deployment} at altitude = { apogee_altitude} and time = {apogee_time}")

# Rocket parameters
# m0 = 50.0 # initial mass [kg]
# mf = 20.0 # final mass [kg]
# burn_time = 20.0 # period of flight w/thrust [s]
# thrust = ve * (m0 - mf)/ burn_time # constant thrust [N]

# # Parachute / Apogee Detection
# if velocity[i-1] > 0 and velocity[i] <= 0:
   # apogee_index = 1
   # apogee_time = t[i]
   # apogee_altitude = altitude[i]
   # parachute_deployment = True
   # print(f"Parachute deployed: {parachute_deployment} at altitude = { apogee_altitude} and time = {apogee_time}")