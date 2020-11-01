function data = ImpactVelocityNotParallelMatlab
    h = waitbar(0, 'Please wait...');

    ThrustTable = xlsread('F15_Thrust','Sheet1');
    mass = 1.0565;
    ThrustAcceleration = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0] / mass);
    ThrustAcceleration = ThrustAcceleration(0:0.01:100);
    waitbar(0, h, 'Loaded Thrust');
    % 
    v = 50:-5:-100;
    numv = length(v);
    s = 1:5:500; % S cannot start at 0, as event will not trigger on first step
    nums = length(s);
    numtotal = nums*numv;
    waitbar(0, h, 'Loaded parameters');

    % [v,s] = meshgrid(v1,s1); %TAKE UP TOO MUCH MEMORY
    % 
    % figure('Visible',true);
    % surface = surf(v, s, NaN(size(v)), 'EdgeColor', 'none');
    % xlabel('Initial Velocity (m/s)');
    % ylabel('Initial Height (m)');
    % colormap(flipud(winter))

    %Q = parallel.pool.DataQueue;c
    %afterEach(Q,@(data) updatePlot(surface,data));

    % Build a waitbar with a cancel button, using appdata to track
    % whether the cancel button has been pressed.

    data = zeros(nums, numv);
    options = odeset('Events',@events,'Refine',4, 'RelTol', 1e-4, 'AbsTol', 1e-6);

    m = 0;
    Q = parallel.pool.DataQueue;
    afterEach(Q, @nUpdateWaitbar);
    waitbar(0, h, 'Initialised');
    
    parfor j = 1:numtotal
        tstart = 0;
        y0 = [s(1 + mod(nums - 1 + j, nums)); v(1 + floor((j - 1)/nums))];

        yeout = [];
        while tstart < 3.45
             [t,y] = ode45(@(t,y) (y(1) > 0 || ThrustAcceleration(1 + round(t * 100)) - 9.80655 > 0) * [y(2); ThrustAcceleration(1 + round(t * 100)) - 9.80655], [tstart 100], y0, options);
             %[t,y] = ode45(@(t,y) ~(y(1) <= 0 && ThrustAcceleration(t) <= g) * [y(2); ThrustAcceleration(t) - g], [tstart tfinal], y0, options);
              %[t,~,~,ye,~] = ode45(@(t,y) ~(y(1) <= 0 && (ThrustAcceleration(1 + round(t / 0.01)) - 9.80655) <= 0) * [y(2); (ThrustAcceleration(1 + round(t / 0.01)) - 9.80655)], [tstart tfinal], y0, options);
    %             [t,y] = ode45(@(t,y) ydot(t, y, ThrustAcceleration), [tstart tfinal], y0, options);
            %yeout = [yeout; y(end)];
            yeout = [yeout; y(end)];
            y0 = [0;0];
            tstart = t(end);
        end
        data(j) = max(abs(yeout));
        send(Q, 1);
    end

%     save 'ImpactVelocity' 'data'
    delete(h)
%     surf(v, s, data, 'EdgeColor', 'none');
%     xlabel('Initial Velocity (m/s)');
%     ylabel('Initial Height (m)');
%     colormap(flipud(winter));

    function nUpdateWaitbar(~)
        m = m + 1;
        waitbar(m/numtotal, h, sprintf('%d / %d', m, numtotal))
%         waitbar(m/numtotal, h);
    end
end

    
% function dydt = ydot(t, y, ThrustAcceleration)
%     a = ThrustAcceleration(1 + round(t / 0.01)) - 9.80655; % MAKE SURE THIS IS ACCURATE
%     dydt = ~(y(1) <= 0 && a <= 0) * [y(2); a]; % IF on the ground and acceleration downwards ignore
% end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(~,y)
    % Locate the time when height passes through zero
    % and stop integration.
    value = y(1) > 0;     % detect height > 0, when height <= 0 gives 0
    isterminal = 1;   % stop the integration
    direction = 0;   % negative direction
end