% close all
% 
ThrustTable = xlsread('F15_Thrust','Sheet1');
Thrust = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0]);
% 
v = 50:-0.01:-100;
numv = length(v);
s = .1:.01:500; % S cannot start at 0, as event will not trigger on first step
nums = length(s);
numtotal = nums*numv;
    
% [v,s] = meshgrid(v1,s1); %TAKE UP TOO MUCH MEMORY
% 
% figure('Visible',true);
% surface = surf(v, s, NaN(size(v)), 'EdgeColor', 'none');
% xlabel('Initial Velocity (m/s)');
% ylabel('Initial Height (m)');
% colormap(flipud(winter))

% Q = parallel.pool.DataQueue;
% afterEach(Q,@(data) updatePlot(surface,data));

% Build a waitbar with a cancel button, using appdata to track
% whether the cancel button has been pressed.
h = waitbar(0, 'Please wait...');

data = zeros(nums, numv);
for j = 1:nums
    for i = 1:numv
        m = nums * (j-1) + i;
        data(j, i) = OneD_case(s(j), v(i), Thrust);
        waitbar(m / numtotal, h, sprintf('%d / %d', m, numtotal))
    end
end
delete(h)

% step = 20;
% partitions = [1:step:numv*nums, numv*nums + 1];
% num_partitions = numel(partitions) - 1;
% 
% clear f;
% f(1:num_partitions) = parallel.FevalFuture;
% 
% % parameterSweep(partitions(1),partitions(2),s,v, Thrust)
% 
% for ii = 1:num_partitions
% %     f(ii) = parfeval(@parameterSweep,1,partitions(ii),partitions(ii+1),s,v, Thrust, Q);
%     f(ii) = parfeval(@parameterSweep,1,partitions(ii),partitions(ii+1),s,v, Thrust);
% end
% 
% updateWaitbarFuture = afterEach(f, @(~) waitbar(sum(strcmp('finished', {f.State}))/numel(f), h, sprintf('%d / %d',sum(strcmp('finished', {f.State})), numel(f))), 1);
% closeWaitbarFuture = afterAll(updateWaitbarFuture, @(h) delete(h), 0);
% 
% data = reshape(fetchOutputs(f), [nums numv]);
save 'ImpactVelocity' 'data'
surface = surf(v, s, data, 'EdgeColor', 'none');
xlabel('Initial Velocity (m/s)');
ylabel('Initial Height (m)');
colormap(flipud(winter))
% 
% % function results = parameterSweep(first,last,s,v, Thrust, Q)
% function results = parameterSweep(first,last,s,v, Thrust)
%     results = zeros(last - first,1);
%     for ii = first:last - 1
% %         result = OneD_case(s(mod(ii, nums)), v( 1 + floor((ii - 1)/nums)), Thrust);
% %         send(Q,[ii,result]);
%         results(ii-first+1) = OneD_case(s(mod(ii, nums)), v( 1 + floor((ii - 1)/nums)), Thrust);
%     end
% end

% function updatePlot(surface,data)
%     surface.ZData(data(1)) = data(2);
%     drawnow('limitrate');
% end