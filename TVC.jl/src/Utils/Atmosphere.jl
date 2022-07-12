using Parameters, Geophysics

export Atmosphere

@with_kw struct Atmosphere @deftype Function
    g = height -> [0; 0; -gravity(0)] #[0; 0; -gravity(height)]
    
    density = Geophysics.density
    speedOfSound = sonicspeed
    kinematicViscocity = kinematic
end