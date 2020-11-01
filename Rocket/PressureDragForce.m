function cd = PressureDragForce(sinphi)
    %{
    persistent a
    
    if isempty(a)
        %%conicalPolyInterpolator
        matrix = [1, 1, 1, 1; 1.3^3, 1.3^2, 1.3, 1; 3, 2, 1, 0; 3*1.3^3, 2*1.3, 1, 0];
        
        cdmach1 = 2.1*sinphi^2 + 0.6019*sinphi;
        y = [sinphi; cdmach1; 4 / (SpecificHeatRatioAir + 1) *(1 - .5*cdmach1); -1.1341 * sinphi];
        a = matrix \ y
    end
    
    input = Mach * FrontalArea / Reference_Area;
    cd = dot(a, [input^3, input^2, input, 1]);
    %}
    cd = .8* sinphi^2;
end