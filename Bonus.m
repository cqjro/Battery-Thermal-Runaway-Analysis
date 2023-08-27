% CHE210 Project
% Group: 7
% Names: Cairo Cristante (1008301348)


% Bonus (+3 marks): Corresponding to Figure 7 [3], 
% show the self-heating rate profiles for lithium cells heated to 
% 80, 100, and 115 °C.

%=========================================================================%
%                              Parameters   
%
%  a -> specific surface area of the sample m^2/g
%  a0 -> reference specific surface area, m^2/g
%  A1 -> frequency factor for the conversion metastable SEI to stable SEI
%        min^-1
%  A2 -> frequency factor for the reaction of intercalated Li with 
%        electrolyte, min^-1
%  aL, aH -> specific surface area of the low and high surface area 
%            materials, respectively, m^2/g
%  C -> specific heat of the ARC sample, J/(g K)
%  d0 -> initial SEI thickness, m
%  dt -> time increment, min
%  dT/dt -> self-heating rate, K/min
%  dxf -> change in the amount of metastable SEI, unitless
%  dxi -> change in the amount of intercalated lithium, unitless
%  E1 -> activation energy for the conversion of metastable SEI to stable 
%        SEI, eV
%  E2 -> activation energy for the reaction of intercalated lithium with 
%        electrolyte, eV
%  h1 -> heat of reaction for the formation of stable SEI from metastable 
%        SEI, J/g
%  h2 -> heat of reaction of intercalated lithium with solvent to 
%        eventually form stable SEI, J/g
%  k1 -> rate constant for the conversion of metastable SEI to stable SEI, 
%        min-1 
%  k2 -> rate constant for the reaction of intercalated lithium with 
%        electrolyte, min-1 
%  kB -> Boltzmann’s constant, eV/K
%  m -> reaction order for the reaction of intercalated lithium with 
%       electrolyte, unitless
%  n -> reaction order for the conversion of metastable SEI to stable SEI, 
%       unitless
%  z -> amount of lithium in the SEI per unit surface area, m-2
%  z0 -> amount of lithium in the SEI per unit surface area, m-2
%=========================================================================%

% Parameter Values
% Note: All parameter values were determed based on the values found in
% Richard et all 1999 II and their experiments.
xf0 = 0.10;
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0 1000]; % larger tpsn than this results in very long run time


% Sets the ODE45 options for this scenario
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);

%% 80ºC Start
T0 = 80;
tspanB = [0:0.1;1000];
y0 = [xf0, xi0(1), T0, z0];
[t1b, y1b] = ode45(@(t, y) scenario3E7b(t, y, aL), tspanB, y0, options);


dydt1b = zeros(length(y1b), 4);
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL);
end

%% 100ºC Start
T0 = 100;
tspanB = [0:0.1:1000];
y0 = [xf0, xi0(1), T0, z0];

[t2b, y2b] = ode45(@(t, y) scenario3E7b(t, y, aL), tspanB, y0, options);

dydt2b = zeros(length(y2b), 4);
for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL);
end


%% 115ºC Start
T0 = 115;
tspanB = [0:0.1:1000];
y0 = [xf0, xi0(1), T0, z0];
[t3b, y3b] = ode45(@(t, y) scenario3E7b(t, y, aL), tspanB, y0, options);

dydt3b = zeros(length(y3b), 4);
for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL);
end
%% Plot
% Self Heating Rate Profile Plot
subplot(2, 2, 1)
semilogy(y1b(:, 3), dydt1b(:, 3))
hold on

semilogy(y2b(:, 3), dydt2b(:, 3))
semilogy(y3b(:, 3), dydt3b(:, 3))
hold off
xlabel('Temperature (ºC)')
ylabel('log dT/dt (ºC/s)')
xlim([80 220])
ylim([10^-2 10^2])
title('Temperature Profile')
legend('7b 80ºC',  '7b 100ºC',  '7b 115ºC')

% Change in the amount of intercalated Li plot
subplot(2, 2, 2)
plot(y1b(:, 3), y1b(:, 2))
hold on
plot(y2b(:, 3), y2b(:, 2))
plot(y3b(:, 3), y3b(:, 2))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of intercalated Li')
xlim([80 220])
title('Change in intercalated Li')
legend( '7b 80ºC',  '7b 100ºC',  '7b 115ºC')

% Change in the amount of metastable SEI plot
subplot(2, 2, 3)
plot(y1b(:, 3), y1b(:, 1))
hold on
plot(y2b(:, 3), y2b(:, 1))
plot(y3b(:, 3), y3b(:, 1))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of Metastable SEI')
xlim([80 220])
title('Change in Metastable SEI')
legend( '7b 80ºC',  '7b 100ºC',  '7b 115ºC')

% Growth of the SEI Layer Plot
subplot(2, 2, 4)
plot(y1b(:, 3), y1b(:, 4))
hold on
plot(y2b(:, 3), y2b(:, 4))
plot(y3b(:, 3), y3b(:, 4))
hold off
xlabel('Temperature (ºC)')
ylabel('Growth of SEI layer')
xlim([80 220])
ylim([0 0.75])
title('Growth of the SEI Layer')
legend( '7b 80ºC',  '7b 100ºC',  '7b 115ºC')


%% Functions

function dydt = scenario3E7b(t, y, a)
    % Parameter Values
    % Note: All parameter values were determed based on the values found in
    % Richard et all 1999 II and their experiments.
    h1C = 150;
    h2C = 325;
    A1 = 1*10^17;
    A2 = 1*10^8;
    E1 = 1.4;
    E2 = 0.8;
    kB = 8.617*10^(-5);
    z0 = 0.15;
    n = 0.5; 
    m = 1; 
    a0 = 1;

% Defines Differential Equations
    dydt = zeros(4,1);
    
    % Solves for xf 
    dydt(1) = -A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n; 
    
    % Solves for xi
    dydt(2) = -A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*exp(-y(4)/z0);
    
    % Solves for T
    dydt(3) = h1C*A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n ...
        + h2C*A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*exp(-y(4)/z0);

    % Solves for z
    dydt(4) = A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*exp(-y(4)/z0);
end

