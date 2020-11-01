numv = 1001;
v = linspace(0, -100,numv);
nums = 101;
s = linspace(0, 100, nums);
numSims = numv * nums;
model = 'ImpactVelocityFunction';

cd 'C:\Users\rag\Documents\GitHub\EPQ\Rocket'
load_system(model)
s0 = 0;
v0 = 0;

for j = 1:numv
    for i = 1:nums
        m = nums * (j-1) + i;
        in(m) = Simulink.SimulationInput(model);
        in(m) = in(m).setVariable('v0', v(j));
        in(m) = in(m).setVariable('s0', s(i));
    end
end

simOut = parsim(in,'ShowSimulationManager','on', 'RunInBackground', 'on', 'UseFastRestart', 'on');

close all

data = zeros(numv, nums);
for j = 1:numv
    for i = 1:nums
            m = nums * (j - 1) + i;
            data(j, i) = max(simOut(m).logsout{1}.Values);
    end
end   

surf(v, s, data')
xlabel('Initial Velocity (m/s)');
ylabel('Initial Height (m)');
colormap(flipud(winter))

h = heatmap(v, s, data, 'Colormap', flipud(gray));
h.Title = 'Impact Velocity';
h.XLabel = 'Initial Velocity (m/s)';
h.YLabel = 'Initial Height (m)';
h.YDisplayData = fliplr(s);
h.ColorLimits = [0 5];