function data = ImpactVelocityParforCompile(v, s, ThrustAcceleration, ThrustStep, tfinal, burnoutTime, M)
    
    numv = length(v);
    nums = length(s);
    numtotal = nums*numv;

    data = zeros(nums, numv);
    options = odeset('Events', @events,'Refine', 4, 'RelTol', 1e-6, 'AbsTol', 1e-9);
    
    parfor (j = 1:numtotal, M)
%     for j = 1:numtotal
        tstart = 0;
        y0 = [s(1 + mod(nums - 1 + j, nums)); v(1 + floor((j - 1)/nums))];

        yeout = [];
        while tstart < burnoutTime
            [t,y] = ode45(@(t,y) (y(1) > 0 || ThrustAcceleration(1 + round(t / ThrustStep)) - 9.80655 > 0) * [y(2); ThrustAcceleration(1 + round(t / ThrustStep)) - 9.80655], [tstart tfinal], y0, options);
            yeout = [yeout; y(end)];
            y0 = [0;0];
            tstart = t(end);
        end
        data(j) = max(abs(yeout));
    end
end

function [value,isterminal,direction] = events(~,y)
    % Locate the time when height passes through zero
    % and stop integration.
    value = double(y(1) > 0);     % detect height > 0, when height <= 0 gives 0
    isterminal = 1;   % stop the integration
    direction = 0;   % negative direction
end