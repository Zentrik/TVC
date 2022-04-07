using Parameters, Geophysics

@with_kw struct atmosphere @deftype Function
    gravity = height -> [0; 0; -gravity(height)]
    
    density = Geophysics.density
    speedOfSound = sonicspeed
    kinematicViscocity = kinematic
end