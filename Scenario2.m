% CHE210 Project
% Group: 7
% Names: Cairo Cristante (1008301348)


% Scenario 2: Estimate the self-heating rate profiles for low-surface area 
% MCMB with different amounts of intercalated Li, showing the self-heating 
% rate profiles calculated based on Eq. 7b, the change in the amount of 
% intercalated Li, the change in the amount of metastable SEI components, 
% and the growth of the SEI layer.

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
T0 = 80;
xf0 = 0.10;
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1; % low surface area MCMB area 1.1 m2/g based on Richard et al I
tspan = 0:0.1:1000;


% Sets the ODE45 options for this scenario
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);

%% 0.0V ODEs solution
y0 = [xf0, xi0(1), T0, z0];
[t1, y1] = ode45(@(t, y) scenario2(t, y), tspan, y0, options);

% 0.0V changes
dydt1 = zeros(length(y1), 4);
for i = 1:length(y1)
   dydt1(i, :) = scenario2(t1(i), y1(i, :));
end

%% 0.089V ODEs solution
y0 = [xf0, xi0(2), T0, z0];
[t2, y2] = ode45(@(t, y) scenario2(t, y), tspan, y0, options);

% 0.089V changes
dydt2 = zeros(length(y2), 4);
for i = 1:length(y2)
   dydt2(i, :) = scenario2(t2(i), y2(i, :));
end

%% 0.127V ODEs solution
tspan = 0:0.1:1500;
y0 = [xf0, xi0(3), T0, z0];
[t3, y3] = ode45(@(t, y) scenario2(t, y), tspan, y0, options);

% 0.127V changes
dydt3 = zeros(length(y3), 4);
for i = 1:length(y3)
   dydt3(i, :) = scenario2(t3(i), y3(i, :));
end

%% Plot
% Self Heating Rate Profile Plot
subplot(2, 2, 1)
semilogy(y1(:, 3), dydt1(:, 3))
hold on
semilogy(y2(:, 3), dydt2(:, 3))
semilogy(y3(:, 3), dydt3(:, 3))
hold off
xlabel('Temperature (ºC)')
ylabel('log dT/dt (ºC/min)')
xlim([80 240])
ylim([10^-2 10^1])
title('Temperature Profile for different intital intecalated Li')
legend('0.0V: x_i = 0.75', '0.089V: x_i = 0.45', '0.127V: x_i = 0.25')

% Change in the amount of intercalated Li plot
subplot(2, 2, 2)
plot(y1(:, 3), y1(:, 2))
hold on
plot(y2(:, 3), y2(:, 2))
plot(y3(:, 3), y3(:, 2))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of intercalated Li')
xlim([40 240])
title('Change in intercalated Li for different intital intecalated Li')
legend('0.0V: x_i = 0.75', '0.089V: x_i = 0.45', '0.127V: x_i = 0.25')

% Change in the amount of metastable SEI plot
subplot(2, 2, 3)
plot(y1(:, 3), y1(:, 1))
hold on
plot(y2(:, 3), y2(:, 1))
plot(y3(:, 3), y3(:, 1))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of Metastable SEI ')
xlim([40 240])
title('Change in Metastable SEI for different intital intecalated Li')
legend('0.0V: x_i = 0.75', '0.089V: x_i = 0.45', '0.127V: x_i = 0.25')

% Growth of the SEI Layer Plot
subplot(2, 2, 4)
plot(y1(:, 3), y1(:, 4))
hold on
plot(y2(:, 3), y2(:, 4))
plot(y3(:, 3), y3(:, 4))
hold off
xlabel('Temperature (ºC)')
ylabel('Growth of SEI layer')
xlim([40 240])
title('Growth of the SEI Layer for different intital intecalated Li')
legend('0.0V: x_i = 0.75', '0.089V: x_i = 0.45', '0.127V: x_i = 0.25')


%% Functions

function dydt = scenario2(t, y)
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
    n = 0.5; %??????
    m = 1; %??????
    a0 = 1;
    a = 1.1;

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






