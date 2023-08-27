% CHE210 Project
% Group: 7
% Names: Cairo Cristante (1008301348)
%        Daniel Moore (1008358457)


% Scenario 4: Provide one feasible alternative (supported by calculations) 
% to mitigate self-heating, “maintaining” the sample temperature at 30 °C 
% (for example with liquid or air cooling).

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
clc

% Parameter Values
% Note: All parameter values were determed based on the values found in
% Richard et all 1999 II and their experiments.
T0 = 80;
xf0 = 0.10;
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1; % low surface area MCMB area 1.1 m2/g based on Richard et al I
h = [600, 1000, 2000];
Tinf = 25;
tspan = 0:0.1:1000; % larger tpsn than this results in very long run time


% Sets the ODE45 options for this scenario
options = odeset('RelTol',1e-5, 'AbsTol', 1e-8, 'NonNegative',[1 2 3 4]);

%% h = 600 J/min m^2
% Solves ODE
y0 = [xf0, xi0(1), T0, z0];
[t1, y1] = ode45(@(t, y) scenario2(t, y, Tinf, h(1)), tspan, y0, options);

%% h = 1000 J/min m^2
% Solves ODE
y0 = [xf0, xi0(1), T0, z0];
[t2, y2] = ode45(@(t, y) scenario2(t, y, Tinf, h(2)), tspan, y0, options);

%% h = 2000 J/min m^2
% Solves ODE
y0 = [xf0, xi0(1), T0, z0];
[t3, y3] = ode45(@(t, y) scenario2(t, y, Tinf, h(3)), tspan, y0, options);

%% Plot
% Self Heating Rate Profile Plot
plot(t1, y1(:, 3))
hold on
plot(t2, y2(:, 3))
plot(t3, y3(:, 3))
hold off
xlabel('Time (min)')
xlim([0, 20])
ylim([15 90])
ylabel('T (ºC)')
title('Temperature Profile')
legend('h = 600 J/min*m^2', 'h = 1000 J/min*m^2', 'h = 2000 J/min*m^2')
%% Functions

function dydt = scenario2(t, y, Tinf, h)
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
    a = 1.1;
    As = (2*pi*12*6+2*pi*12^2)/10^6; % m^2 surface area of battery
    % The battery is assumed to be a 25mm diameter, 6 mm thick
    mass = 3/1000; % kg assume 1 gram battery
    Cp = 830; % J/kg*K specific heat capacity of battery
    b = h*As/(mass*Cp); % 1/tau
    

% Defines Differential Equations
    dydt = zeros(4,1);
    
    % Solves for xf 
    dydt(1) = -A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n; 
    
    % Solves for xi
    dydt(2) = -A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*exp(-y(4)/z0);
    
    % Solves for T
    dydt(3) = h1C*A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n ...
        + h2C*A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*exp(-y(4)/z0) ...
        - b*(y(3) - Tinf); 
    % follows lumped system cooling through b(T(t) - Tinf)

    % Solves for z
    dydt(4) = A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*exp(-y(4)/z0);
end






