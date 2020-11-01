function impactVelocity = OneD_case(s0, v0, Thrust)
    tstart = 0;
    tfinal = 30;
    % plot(-2:0.001:6, Thrust(-2:0.001:6))
    y0 = [s0; v0];
    options = odeset('Events',@events,'Refine',4, 'RelTol', 1e-4, 'AbsTol', 1e-6);

%        tout = tstart;
%        yout = y0.';
    % teout = [];
    yeout = [];
    % ieout = [];
    
    while tstart < 3.45
        % Solve until the first terminal event.
%         [t,y,te,ye,ie] = ode45(@(t,y) ydot(t, y, Thrust), [tstart tfinal], y0, options);
        [t,~,~,ye,~] = ode45(@ydot, [tstart tfinal], y0, options);
        % Accumulate output.  This could be passed out as output arguments.
        %nt = length(t);
        %ye_t = [ye_t ye(2, 2)];
%            tout = [tout; t(2:end)];
%            yout = [yout; y(2:end,:)];
    %     teout = [teout; te];          % Events at tstart are never reported.
         yeout = [yeout; ye];
    %     ieout = [ieout; ie];

        % Set the new initial conditions, with .9 attenuation.
        y0 = [0;0];
        tstart = t(end);
    end
    impactVelocity = max(abs(yeout(:,2)));
    
    function dydt = ydot(t, y)
        a = Thrust(t)/1.0565 - 9.80655;
        if y(1) <= 0 && a <= 0
            dydt = [0;0];
        else
            dydt = [y(2); a];
        end
    end
end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y)
    % Locate the time when height passes through zero
    % and stop integration.
    value = double(y(1) > 0);     % detect height > 0, when height <= 0 gives 0
    isterminal = 1;   % stop the integration
    direction = 0;   % negative direction
end