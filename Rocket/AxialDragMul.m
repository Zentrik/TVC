function mul = AxialDragMul(AoA)
    persistent sn  p2  a a2
    
    if isempty(sn)
        sn = deg2rad(17);
        p2 = deg2rad(90);

        %% AxialPoly1
        matrix = [0, 0, 0, 1; sn^3, sn^2, sn, 1; 0 ,0, 1, 0; 3*sn^2, 2*sn, 1, 0];
        a = matrix \ [1;1.3;0;0];

        %% AxialPoly2
        matrix2 = [sn^4, sn^3, sn^2, sn, 1 ;p2^4, p2^3, p2^2, p2, 1 ;4*sn^3, 3*sn^2, 2*sn, 1, 0 ;4*p2^3, 3*p2^2, 2*p2, 1, 0 ;12*p2^2, 6*p2, 2, 0, 0];
        a2 = matrix2 \ [1.3;0;0;0;0];
    end
    
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
end