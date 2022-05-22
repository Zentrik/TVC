using LinearAlgebra
using ..Utils

export CalculateAero

function CalculateAero(x, p, t)
    veh = p.veh
    atmos = p.atmos
    wind = p.wind

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]
    
    Airspeed = rotate(conjugate(quat), v + wind)

    lengtha = norm(Airspeed[1:2])

    if lengtha ≈ 0 
        AeroForce = zeros(3)
        AeroTorque = zeros(3)
    else
        DynamicPressure = atmos.density(height) * norm(Airspeed)^2 / 2
        AoA = acos(Airspeed[3] / norm(Airspeed))
        Mach = norm(Airspeed) / atmos.speedOfSound(height)
        
        cosa = Airspeed[1] / lengtha;
        sina = Airspeed[2] / lengtha;

        ThetaRotation = [cosa -sina 0; sina cosa 0; 0 0 1]

        ω_local = ThetaRotation \ ω # ω in wind coordinates?

        (CN, Cpitch) = PitchNormalCD(veh, Mach, AoA, ω_local[2], Airspeed, t)
        Caxial = CalculateCD(veh, atmos, Mach, AoA, Airspeed, height)
        Cyaw = 0
        Croll = 0
        Cside = 0
            
        AeroForce = ThetaRotation * DynamicPressure * veh.Reference_Area * -[CN; Cside; Caxial] # convert aero forces in wind frame to body frame?
        AeroTorque = ThetaRotation * (DynamicPressure * veh.Reference_Area * veh.Reference_Diameter * [-Cyaw; Cpitch; Croll])
    end

    return (force=AeroForce, torque=AeroTorque)
end

function PitchNormalCD(rocket, Mach, AoA, pitchrate, Airspeed, t)
    if Mach < 0.05 && AoA > pi/4 
        cnmul = (Mach/0.05)^2;
    else
        cnmul = 1.0
    end

    NoseCone_CN = sin(AoA) * ( sin(AoA) * cnmul * rocket.K * rocket.NoseCone_PlanformArea / rocket.Reference_Area + 2 );
    BodyTube_CN = sin(AoA) * ( sin(AoA) * cnmul * rocket.K * rocket.BodyTube_PlanformArea / rocket.Reference_Area );
    NormalForceCD = BodyTube_CN + NoseCone_CN;

    CNCPxL = ( NoseCone_CN * rocket.NoseCone_CP + BodyTube_CN * rocket.BodyTube_CP) / rocket.Reference_Diameter;

    cacheDiameter = (rocket.NoseCone_PlanformArea + rocket.BodyTube_PlanformArea) / (rocket.NoseCone_Length + rocket.BodyTube_Length);

    mul = 3 * (rocket.PitchCenterX^4 + (rocket.NoseCone_Length + rocket.BodyTube_Length - rocket.PitchCenterX)^4 ) * (.275 * cacheDiameter / (rocket.Reference_Area * rocket.Reference_Diameter));
    PitchDampingMomentCD = min( ( mul * ( pitchrate / norm(Airspeed) )^2 ) , CNCPxL ) * sign(pitchrate);

    CNCmxL = NormalForceCD * rocket.CG(t) / rocket.Reference_Diameter;
    PitchMomentCD = CNCPxL - CNCmxL - PitchDampingMomentCD;
    return (NormalForceCD, PitchMomentCD)
end

function CalculateCD(rocket, atmos, Mach, AoA, Airspeed, height)
    Cbase = .12 + .13*Mach^2;
    CFrictionDrag = CalculateFrictionDrag(rocket, atmos, Airspeed, Mach, height)
    Cpressure = PressureDragForce(rocket);
    
    Cdrag = Cbase + CFrictionDrag + Cpressure;
    Caxial = AxialDragMul(AoA) * Cdrag;
    return Caxial
end

function CalculateFrictionDrag(rocket, atmos, Airspeed, Mach, height)
    Reynolds = norm(Airspeed) * (rocket.NoseCone_Length + rocket.BodyTube_Length) / atmos.kinematicViscocity(height);
    
    Cf = (1.5 * log(Reynolds) - 5.6)^-2;
    c1 = 1 - 0.1 * Mach^2;
    Cfc1 = Cf * c1;
    
    RoughnessCorrection = c1;
    NoseCone_RoughnessLimited = .032 * (rocket.SurfaceRoughness / rocket.NoseCone_Length)^0.2 * RoughnessCorrection;
    BodyTube_RoughnessLimited = .032 * (rocket.SurfaceRoughness / rocket.BodyTube_Length)^0.2 * RoughnessCorrection;
    
    NoseCone_Cf = max(Cfc1, NoseCone_RoughnessLimited);
    BodyTube_Cf = max(Cfc1, BodyTube_RoughnessLimited);
    
    BodyFriction = BodyTube_Cf * rocket.BodyTube_WetArea + NoseCone_Cf * rocket.NoseCone_WetArea;
    
    fB = (rocket.NoseCone_Length + rocket.BodyTube_Length + 0.0001) / rocket.Reference_Radius;
    correction = 1 + 1 / (2 * fB);
    
    FrictionCd = BodyFriction * correction / rocket.Reference_Area
    return FrictionCd
end

function PressureDragForce(rocket)
   return 0.8 * rocket.sinphi^2;
end

function AxialDragMul(AoA)    
    sn = deg2rad(17);
    p2 = deg2rad(90);

    # AxialPoly1
    matrix = [0 0 0 1; sn^3 sn^2 sn 1; 0 0 1 0; 3*sn^2 2*sn 1 0];
    a = matrix \ [1; 1.3; 0; 0];

    # AxialPoly2
    matrix2 = [sn^4 sn^3 sn^2 sn 1; p2^4 p2^3 p2^2 p2 1; 4*sn^3 3*sn^2 2*sn 1 0; 4*p2^3 3*p2^2 2*p2 1 0; 12*p2^2 6*p2 2 0 0]
    a2 = matrix2 \ [1.3; 0; 0; 0; 0];
    
    AoAcalc = AoA;
    
    if AoAcalc > p2
        AoAcalc = pi - AoAcalc;
    end

    if AoAcalc < sn
        mul = dot(a, [AoAcalc^3, AoAcalc^2, AoAcalc, 1]);
    else 
        mul = dot(a2, [AoAcalc^4, AoAcalc^3, AoAcalc^2, AoAcalc, 1]);
    end

    if AoA >= p2
        mul = -mul;
    end

    return mul
end