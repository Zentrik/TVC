import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import quaternion

thrustcsv = np.genfromtxt("OpenRocket/OpenRocket Thrust F10.csv", delimiter=',')
CenterOfMasscsv = np.genfromtxt("Rocket/mass,cg,I.csv", delimiter=',')

K = 1.1
Mach = 0

Reference_Diameter = 7.62e-2
Reference_Radius = Reference_Diameter / 2
Reference_Area = Reference_Radius^2 * np.pi

NoseCone_Length = 15e-2
BodyTube_Length = 80e-2
NoseCone_WetArea = np.pi * Reference_Radius * NoseCone_Length          # np.pi*r*l 
BodyTube_WetArea = np.pi * Reference_Radius * 2 * BodyTube_Length      # np.pi*d*l
NoseCone_PlanformArea = Reference_Diameter * NoseCone_Length / 2
BodyTube_PlanformArea = Reference_Diameter * BodyTube_Length

SurfaceRoughness = 20e-6
    
sinphi = Reference_Radius / np.hypot(Reference_Radius, NoseCone_Length)
BodyTube_CP = 0.55
NoseCone_CP = 0.1
pitchCenterX = 0


airSpeed = np.zeros(3)
orientation = np.quaternion(1,0,0,0)

timestep = 0.0001
timerange = np.arange(0.0, 6.0, timestep)

def calculateAcceleration(thrust, CenterOfMassX, COT_to_COM, mass, Ixx, orientation, globalAngularVelocity):
    airSpeedLength = np.linalg.norm(airSpeed)
    AoA = np.arccos(airSpeed[2] / airSpeedLength)

    if Mach < 0.05 and AoA > np.pi/4:
        cnmul = (Mach/0.05)^2
    else:
        cnmul = 1

    NoseCone_CN = np.sin(AoA) * ( np.sin(AoA) * cnmul * K * NoseCone_PlanformArea / Reference_Area + 2 )
    BodyTube_CN = np.sin(AoA) * ( np.sin(AoA) * cnmul * K * BodyTube_PlanformArea / Reference_Area )
    NormalForceCD = BodyTube_CN + NoseCone_CN

    CNCPxL = ( NoseCone_CN * NoseCone_CP + BodyTube_CN * BodyTube_CP) / Reference_Diameter

    cacheDiameter = (NoseCone_PlanformArea + BodyTube_PlanformArea) / (NoseCone_Length + BodyTube_Length)

    mul = 3 * (pitchCenterX^4 + (NoseCone_Length + BodyTube_Length - pitchCenterX)^4 ) * (.275 * cacheDiameter / (Reference_Area * Reference_Diameter))
    pitchDampingMomentCD = min( ( mul * ( pitchrate / airSpeedLength )^2 ) , CNCPxL )   * sign(pitchrate)

    CNCmxL = NormalForceCD * CenterOfMassX / Reference_Diameter
    pitchMomentCD = CNCPxL - CNCmxL - pitchDampingMomentCD

    acceleration = thrust / mass
    if Xe[2] > 0:
        angularAcceleration = np.array([0, pitchMomentCD, 0]) / Ixx
    else:
        angularAcceleration = np.zeros(3)

    Ae = quaternion.rotate_vectors(orientation, acceleration)
    globalAngularAcceleration = quaternion.rotate_vectors(orientation, angularAcceleration)

    return Ae, globalAngularAcceleration

for t in timerange:
    thrust = np.interp(t, thrustcsv[:,0], thrustcsv[:,1])
    CenterOfMassX = np.interp(t, CenterOfMasscsv[:,0], CenterOfMasscsv[:,4]) / 100
    COT_to_COM = .95 - CenterOfMassX
    mass = np.interp(t, CenterOfMasscsv[:,0], CenterOfMasscsv[:,1]) / 1000
    Ixx = np.interp(t, CenterOfMasscsv[:,0], CenterOfMasscsv[:,2])

    k1Ae, k1gAA = calculateAcceleration(thrust, CenterOfMassX, COT_to_COM, mass, Ixx, orientation, globalAngularVelocity)


    # Second position, k2 = f(t + h/2, y + k1*h/2)

    k2Xe = Xe + Ve * timestep / 2
    k2Ve = Ve + k1Ae * timestep / 2
    k2quat =  quat + timestep/2 * .5 * np.quaternion(0, gW[0], gW[1], gW[2]) *  quat
    k2gW = gW + k1gAA * timestep / 2
    
    k2Ae, k2gAA = calculateAcceleration(thrust, CenterOfMassX, COT_to_COM, mass, Ixx, orientation, globalAngularVelocity)
    

    # Third position, k3 = f(t + h/2, y + k2*h/2)

    k3Xe = k2Xe + k2Ve * timestep / 2
    k3Ve = k2Ve + k2Ae * timestep / 2
    k3quat =  k2quat + timestep/2 * .5 * np.quaternion(0, k2gW[0], k2gW[1], k2gW[2]) * k2quat
    k3gW = K2gW + k2gAA * timestep / 2
    
    k3Ae, k3gAA = calculateAcceleration(thrust, CenterOfMassX, COT_to_COM, mass, Ixx, orientation, globalAngularVelocity)
    

    # Fourth position, k4 = f(t + h, y + k3*h)
    
    k4Xe = k1Xe + k3Ve * timestep
    k4Ve = k1Ve + k3Ae * timestep
    k4quat =  k1quat + timestep * .5 * np.quaternion(0, k3gW[0], k3gW[1], k3gW[2]) * k3quat
    k4gW = K1gW + k3gAA * timestep
    
    k4Ae, k4gAA = calculateAcceleration(thrust, CenterOfMassX, COT_to_COM, mass, Ixx, orientation, globalAngularVelocity)
    

    # Sum all together,  y(n+1) = y(n) + h*(k1 + 2*k2 + 2*k3 + k4)/6
    Ve = Ve + timestep/ 6 * (k1Ae + 2*k2Ae +2*k3Ae + k4Ae)
    Xe = Xe + timestep/ 6 * (k1Ve + 2*k2Ve +2*k3Ve + k4Ve)
    gW = gW + timestep / 6 * (k1gAA + 2*k2gAA +2*k3gAA + k4gAA)
    deltaTheta = timestep / 6 * (k1gW + 2*k2gW +2*k3gW + k4gW)
    
    quat = quat + np.quaternion(0, deltaTheta[0], deltaTheta[1], deltaTheta[2]) * quat