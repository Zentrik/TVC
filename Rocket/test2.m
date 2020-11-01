gridSize = 40;
sigma = linspace(5, 45, gridSize);
numrho = 50;
rho = linspace(-50, -100, numrho);
beta = 8/3;

[rho,sigma] = meshgrid(rho,sigma);

figure('Visible',true);
surface = surf(rho,sigma,NaN([gridSize numrho]));
xlabel('\rho','Interpreter','Tex')
ylabel('\sigma','Interpreter','Tex')

parpool;

Q = parallel.pool.DataQueue;
afterEach(Q,@(data) updatePlot(surface,data));

step = 100;
partitions = [1:step:gridSize * numrho, gridSize * numrho + 1];

f(1:numel(partitions) - 1) = parallel.FevalFuture;

for ii = 1:numel(partitions) - 1
    f(ii) = parfeval(@parameterSweep,1,partitions(ii),partitions(ii+1),sigma,rho,beta,Q);
end

function results = parameterSweep(first,last,sigma,rho,beta,Q)
    results = zeros(last-first,1);
    for ii = first:last-1
        lorenzSystem = @(t,a) [sigma(ii)*(a(2) - a(1)); a(1)*(rho(ii) - a(3)) - a(2); a(1)*a(2) - beta*a(3)];
        [t,a] = ode45(lorenzSystem,[0 100],[1 1 1]);
        result = rocket(sigma(ii), rho(ii));
        send(Q,[ii,result]);
        results(ii-first+1) = result;
    end
end

function updatePlot(surface,data)
    surface.ZData(data(1)) = data(2);
    drawnow('limitrate');
end

