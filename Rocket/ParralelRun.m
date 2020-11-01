numv = 12;
v = linspace(-3,-8,numv);
nums = 6;
s = linspace(5,7, nums);
numSims = numv * nums;
model = 'roll_control';

cd 'C:\Users\rag\Documents\GitHub\EPQ\Rocket'
load_system(model)
s0 = 0;
v0 = 0;
Simulink.BlockDiagram.buildRapidAcceleratorTarget(model);

for j = 1:numv
    for i = 1:nums
        m = nums * (j-1) + i;
        in(m) = Simulink.SimulationInput(model);
        in(m) = in(m).setVariable('v0',v(j));
        in(m) = in(m).setVariable('s0',s(i));
    end
end

in = in.setModelParameter('SimulationMode', 'rapid-accelerator');
in = in.setModelParameter('RapidAcceleratorUpToDateCheck', 'off');
simOut = parsim(in,'ShowSimulationManager','on', 'TransferBaseWorkspaceVariables', 'on');

close all

rows = nums;
for j = 1:rows
    subplot(rows, 2, 2*j - 1)
    legend_labels = cell(1,numv);
    for m = 1:numv
            i = nums * (m - 1) + j;
            out = simOut(i);
            ts = out.tmp_raccel_logsout{8}.Values;
            data = ts.data(:,3);
            plot(ts.time, data);
            legend_labels{m} = ['Run ' num2str(i)];
            hold all
    end
    xlabel('Time (s)');
    ylabel('Height (m)');
    legend(legend_labels)

    subplot(rows, 2, 2*j)

    for m = 1:numv
            i = nums * (m - 1) + j;
            out = simOut(i);
            ts = out.tmp_raccel_logsout{7}.Values;
            data = ts.data(:,3);
            plot(ts.time, data);
            hold all
    end
    xlabel('Time (s)');
    ylabel('Z Velocity (m)');
    legend(legend_labels)
end