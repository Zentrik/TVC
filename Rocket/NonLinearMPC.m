%% create MPC controller object with sample time
nx = 12;
ny = 12;
nlobj = nlmpc(nx, ny, 'MV', [1 2 3], 'MD', [4 5 6 7 8]);

nlobj.Model.StateFcn = "RocketDerivative";
%% specify Time step
nlobj.Ts = 0.1;
%% specify prediction horizon
nlobj.PredictionHorizon = 200;
%% specify control horizon
nlobj.ControlHorizon = 200;

nlobj.Optimization.CustomCostFcn = @(X,U,e,data) Ts*sum(sum(U(1:p,:)));
nlobj.Optimization.ReplaceStandardCost = true;
nlobj.Optimization.CustomEqConFcn = @(X,U,data) X(end,:)';

%%
x0 = zeros(1, 12);
u0 = [0 0 0];
md0 = [10.6, 0.0826975856, 0.0826975856, 2.4778e-04, 0.401242753755115];
validateFcns(nlobj,x0,u0,md0);