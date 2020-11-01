% data = matfile('ImpactVelocity');
% sz = 1:0.01:500;
% vz = 50:-0.01:-100;
% % % z = data.data;
% % vz = v;
% % sz = s;
% fun = griddedInterpolant({sz, fliplr(vz)}, fliplr(data.data));
% % surf(vz, sz, fun({sz, vz}), 'EdgeColor', 'none');
% % s0 = 30;
% % v0 = -20;
% % gamma = 0:0.001:(-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2);
% % plot(gamma, fun(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma));
% % [x,fval,exitflag] = fminbnd(@(gamma) fun(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma) , 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2));
%  
% numv = length(vz);
% nums = length(sz);
% numtotal = nums*numv;

% gammadata = zeros(nums, numv);
% save OptimalGamma gamma -v7.3
% file = 'OptimalGamma.mat';
% matObj = matfile(file,'Writable',true);
% % Q = parallel.pool.DataQueue;
% % afterEach(Q, @sa);
% % 
% % parfor (h = 1:numtotal, 4)
% % % for h = 1:numtotal
% %     s0 = sz(1 + mod(nums - 1 + h, nums));
% %     v0 = vz(1 + floor((h - 1)/nums));
% % 
% %     og = fminbnd(@(gamma) fun(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma) , 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2));
% %     send(Q, [og, h]);Op
% % end
% % 
% % function sa(a)
% %      OptimalGamma.gamma(a(2)) = a(1);
% % end

% parfor (h = 1:numtotal, 4)
%     s0 = sz(1 + mod(nums - 1 + h, nums));
%     v0 = vz(1 + floor((h - 1)/nums));
% 
%     gammadata(h) = fminbnd(@(gamma) fun(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma) , 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2));
% end

% data = matfile('ImpactVelocity.mat');
% Datasize = size(data, 'data');

% s = 1:0.01:500;
% v = 50:-0.01:-100;
s = 1:0.1:150;
v = 50:-0.1:-50;

% decimated = data.data(1:10:Datasize(1), 1:10:Datasize(2));
% fun = griddedInterpolant({s, fliplr(v)}, fliplr(data.data), 'nearest', 'nearest');

% [smesh, vmesh] = meshgrid(sz, vz);
% fun = fit([smesh(:), vmesh(:)], decimated(:), 'lowess');
%  
% surf(vz, sz, fun({sz, vz}));
% shading interp
% fun = griddedInterpolant({s, fliplr(v)}, fliplr(data.data));

% ThrustTable = xlsread('F15_Thrust','Sheet1');
% mass = 1.0565;
% Thrust = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0]);
% ThrustAccelerationInterpolant = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0] / mass);
% 
% ThrustStep = 1e-3;
% tfinal = 30;
% ThrustAcceleration = ThrustAccelerationInterpolant(0:ThrustStep:tfinal);
% burnoutTime = 3.45;

% gammadata = GenerateData(v, s, ThrustAcceleration, ThrustStep, tfinal, burnoutTime, 4);
% save OptimalGammaDecimated gammadata -v7.3
% surf(v, s, gammadata)
% shading interp

% vz = 0:-0.1:-50;
% sz = 1:0.1:10;
% [~,idxs] = histc(sz,s);
% [vs,s_idx] = sort(v);
% [~,tmp] = histc(vz,vs);
% idxv = s_idx(tmp);
% surf(vz, sz, gammadata(idxs, idxv))
% shading interp

% s0 = 4.2; 
% v0 = -23.4;
% gamma = 0:10e-3:(-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2) - 10e-10;
% x = zeros(1, length(gamma));
% odeoptions = odeset('Events', @events,'Refine', 4, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% for i = 1:length(gamma)
%    x(i) = simulate(s0 + v0.*gamma(i) - 9.80655.*gamma(i).^2/2, v0 - 9.8055.*gamma(i), tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions);
% end
% plot(gamma, x)
% plot(gamma, fun(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma))

% gamma = 0;
% i = 1;
% simulate(s0 + v0.*gamma(i) - 9.80655.*gamma(i).^2/2, v0 - 9.8055.*gamma(i), tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions);

load OptimalGammaDecimated
% [snn, vnn] = ndgrid(s, v);
% inputNN = [reshape(snn, 1, []); reshape(vnn, 1, [])];
% outputNN = reshape(gammadata, 1, []); 
% estimatedgammadata = reshape(gammadataNN(inputNN), size(snn));
% surf(vz, sz, gammadata(idxs, idxv), 'FaceColor','g', 'FaceAlpha',0.5)
% hold on
% surf(vz, sz, estimatedgammadata(idxs, idxv), 'FaceColor','r', 'FaceAlpha',0.5)
% hold off

delayforignitor = 0.1; %time ignitor takes to ignite
samplerate = 0; %time between samples

bound = delayforignitor; % want to ignite when optimal time = delayforigntitor. To deal with discretisation assuming sawtooth shape, you want to ignite now if optimaltime-delayforignitor <= optimaltime
% Ignored discretisation as then position and velocity at next time step
% would have to be estimated to get optimal time at next time step.
% gammadata( gammadata >= bound ) = NaN;
% gammadata( gammadata < bound ) = 0;
% sum(gammadata(:)==0) % number of rows/columns found.
[row,col] = find(gammadata < bound);
whentoignite = [s(row); v(col)]; %indices are coded by s*10 - 9 and v*-10 + 501
% surf(v, s, gammadata)
% shading interp

state = round([10; -24], 1); %current state rounded to 1dp
tol = 5e-2;
if sum(sum(abs(whentoignite - state), 1) < tol) > 0 % take the absolute of subtracted matrix and sum columwise, find number of elements less than tol, if greater than 0 ignite
    1 % ignite
end

function gammadata = GenerateData(v, s, ThrustAcceleration, ThrustStep, tfinal, burnoutTime, M)
    h = waitbar(0, 'Please wait...');

    numv = length(v);
    nums = length(s);
    numtotal = nums*numv;

    gammadata = zeros(nums, numv);

    m = 0;
    Q = parallel.pool.DataQueue;
    afterEach(Q, @nUpdateWaitbar);
    waitbar(0, h, 'Initialised');
    
    options = optimoptions('patternsearch', 'Display', 'off', 'Cache', 'on');
    odeoptions = odeset('Events', @events,'Refine', 4, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    
    waitbar(0, h, 'Starting Parallel Pool');
    parfor (l = 1:numtotal, M)
        s0 = s(1 + mod(nums - 1 + l, nums));
        v0 = v(1 + floor((l - 1)/nums));

%         gammadata(l) = fminbnd(@(gamma) Interpolant(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma) , 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2));
%         simulannealbnd(@(gamma) Interpolant(s0 + v0.*gamma -
%         9.80655.*gamma.^2/2, v0 - 9.8055.*gamma), (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2), 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2));
%         while gammadata(l) == 0 %sometimes 0 is returned randomly when it true answer is not 0
%              gammadata(l) = particleswarm(@(gamma) Interpolant(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma), 1, 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2), options); %faster than annealing, better than fminbnd
%         end
%          simulate(s0, v0, tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions)
%       simulannealbnd(@(gamma) simulate(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma, tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions), (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2), 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2) - eps);
        gammadata(l) = patternsearch(@(gamma) simulate(s0 + v0.*gamma - 9.80655.*gamma.^2/2, v0 - 9.8055.*gamma, tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions), (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2), [], [], [], [], 0, (-v0 - (v0^2-4*s0*-9.80655/2)^0.5)/(2*-9.80655/2) - 10e-10, [], options);
        send(Q, 1);
    end
    
    delete(h)
    function nUpdateWaitbar(~)
        m = m + 1;
        waitbar(m/numtotal, h, sprintf('%d / %d', m, numtotal))
    end  
end   

function maximpact = simulate(s0, v0, tfinal, burnoutTime, ThrustAcceleration, ThrustStep, odeoptions)
    tstart = 0;
    y0 = [s0; v0];

    yeout = [];
%     tout = tstart;
%     yout = y0.';
    while tstart < burnoutTime
        [t,y] = ode45(@(t,y) (y(1) > 0 || ThrustAcceleration(1 + round(t / ThrustStep)) - 9.80655 > 0) * [y(2); ThrustAcceleration(1 + round(t / ThrustStep)) - 9.80655], [tstart tfinal], y0, odeoptions);
        yeout = [yeout; y(end)];
        y0 = [0;0];
%         tout = [tout; t(2:end)];
%         yout = [yout; y(2:end,:)];
        tstart = t(end);
%             if tstart >= tfinal
%                 [s0, v0]
%             end
    end
    maximpact = max(abs(yeout));
%     plot(tout, yout)
end

function [value,isterminal,direction] = events(~,y)
    % Locate the time when height passes through zero
    % and stop integration.
    if y(1) > 0 %&& t ~= tstart % detect height > 0, when height <= 0 gives 0 if not on first timestep
        value = 1;
    else
        value = 0;
    end
%             value = double(y(1) > 0);  
    isterminal = 1;   % stop the integration
    direction = 0;   % negative direction
end