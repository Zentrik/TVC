325.059163748137/mean(masscgI(:,2)) - .5 * g* 7.13^2;
trapz(cumtrapz(ThrustTable(:, 1), ThrustTable(:, 2)))/mean(masscgI(:,2)) - .5 * 7.13^2 * g;

Impulse = 76.3253;
IntegralofImpulse = 325.059163748137;
burntime = 7.13;

Impulse = 74.3897504622088;
IntegralofImpulse = 308.858714637416;
burntime = 7.02;

Impulse = 16.7391811190443;
IntegralofImpulse = 15.5163001126130;
burntime = 1.65;

%% F15
Impulse = 49.2278625190001;
IntegralofImpulse = 87.3413889289923;
burntime = 3.45;

%%Flight Data LANDS .5s too early when using drag
Impulse = 49.1715140996190;
IntegralofImpulse = 87.2121878193817;
DragImpulse = 0.247389841037678;
DragIntegralofImpulse = 0.671316809574649;

DragImpulse = 0.271187253395773;
DragIntegralofImpulse = 0.747767832381193;

DragImpulse = 0.15;
DragIntegralofImpulse = 0.3;

%% Smaller step size
Impulse = 49.2157275000001;
IntegralofImpulse = 87.3494947694980;

g = 9.80655;
mass = mean(masscgI(:,2));
v0 = g*burntime - Impulse/mass - DragImpulse/mass
s0 = -v0 * burntime + .5 * g * burntime^2 - IntegralofImpulse/mass - DragIntegralofImpulse/mass