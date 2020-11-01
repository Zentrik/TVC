num = 100;

model = 'roll_control';

cd 'C:\Users\rag\Documents\GitHub\EPQ\Rocket'
load_system(model)

% plot(simOut(1).tmp_raccel_logsout{64}.Values)
% hold on
% plot(simOut(2).tmp_raccel_logsout{64}.Values)
clear simOut

% K = zeros(3, 9);

warning('off','all')
for i = 1:num
    simOut(i) = sim(model,'SimulationMode','accelerator');
end
warning('on','all')

data = [];
int = 0;
b = 0;
for i = 1:num   
    simdata = reshape(simOut(i).tmp_raccel_logsout{2}.Values.Data, [4 simOut(i).tmp_raccel_logsout{2}.Values.Length]);
    simdata(4, :) = rad2deg(simdata(4, :));
    data = [data simdata];
    int = int + trapz(simOut(i).tmp_raccel_logsout{2}.Values.Time, [simdata(1:3, :); abs(simdata(4,:))], 2);
end

mean = int/(3.45 * num);

for i = 1:num
    simdata = reshape(simOut(i).tmp_raccel_logsout{2}.Values.Data, [4 simOut(i).tmp_raccel_logsout{2}.Values.Length]);
    simdata(4, :) = rad2deg(simdata(4, :));
    b = b + trapz(simOut(i).tmp_raccel_logsout{2}.Values.Time, ([simdata(1:3, :); abs(simdata(4,:))] - mean).^2, 2);
end    

variance = b/(3.45 * num);
maximum = [max(abs(data(1:3, :).')) max(abs(data(4, :).'))];

mean.'
variance.'
maximum