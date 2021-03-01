function [y1] = gammadataNN(x1)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 14-Jul-2020 23:08:19.
%
% [y1] = myNeuralNetworkFunction(x1) takes these arguments:
%   x = 2xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [1;-50];
x1_step1.gain = [0.0134228187919463;0.02];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.1665979669494142001;-1.4356811623133816092;0.95782137863719551962;0.64586843994520937162;-1.3083950902043395281;-0.66197983503488644352;-2.0638722866789738219;2.0742340940744377065;5.8528519390417894641;2.0630548747443593349];
IW1_1 = [-1.6632088254344321587 -3.4808636988325334372;0.30968638091687500369 1.1141230434669797678;-2.4961822442387942012 2.5395072913000760195;1.0654971050803523358 0.98283329349031489652;0.25068982251137894579 1.0074700479176561529;-1.0372490490987511524 -1.0116834566633412518;-2.4609741851915880595 -2.7768578252068922829;2.1255620030751916083 2.8475143573413630449;5.3322528350426470212 11.888906698318461252;2.3036079295275206924 2.8019312010846362249];

% Layer 2
b2 = 1.3343912881691786243;
LW2_1 = [0.081126712739063980284 -6.7358964469498996408 -0.0046173080774855384742 -3.2243566848204885389 9.1071598051327029566 -3.537540753035568919 3.2817196707966824754 -3.8420472387685324911 0.025604824799680989816 6.9187286751112404559];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.195339075597689;
y1_step1.xoffset = 1.12789777517719e-11;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end