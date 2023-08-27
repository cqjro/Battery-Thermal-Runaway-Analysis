% CHE210 Project
% Group: 7
% Prepared by: Cairo Cristante (1008301348) 
%              Anastasia Dimov (1006880038) 
%              Miranda Su (1007039973)

% Scenario 3: Corresponding to Figure 5 [3], 
% show the difference between the low- and high- surface-area MCMB 
% calculated using both Eq. 7a and 7b. Comment on the advantages and 
% disadvantages of selecting a low- vs high- surface-area electrode.

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
T0 = 100;
xf0 = 0.10;
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = 1:0.1:1000;

% Sets the ODE options for this scenario
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);


%% Low Surface Area
y0 = [xf0, xi0(1), T0, z0]; % defines inital conditions

% Solves the ODEs in the A and B version of functions
[t1a, y1a] = ode45(@(t, y) scenario3E7a(t, y, aL), tspan, y0, options);
[t1b, y1b] = ode45(@(t, y) scenario3E7b(t, y, aL), tspan, y0, options);

dydt1a = zeros(length(y1a), 4); % intializes accumulators
dydt1b = zeros(length(y1b), 4);
for i = 1:length(y1a) % evaluates the change using solved points in ODE func
   dydt1a(i, :) = scenario3E7a(t1a(i), y1a(i, :), aL);
end
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL);
end

%% High Surface Area
xf0 = 0.16; % From Richard et all II for high surface area
y0 = [xf0, xi0(1), T0, z0];
[t2a, y2a] = ode45(@(t, y) scenario3E7a(t, y, aH), tspan, y0, options);
[t2b, y2b] = ode45(@(t, y) scenario3E7b(t, y, aH), tspan, y0, options);

dydt2a = zeros(length(y2a), 4);
dydt2b = zeros(length(y2b), 4);
for i = 1:length(y2a)
   dydt2a(i, :) = scenario3E7a(t2a(i), y2a(i, :), aH);
end
for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aH);
end


%% Plot
%For changing color of lines
strA = '#EDB120';
colorA = sscanf(strA(2:end),'%2x%2x%2x',[1 3])/255;

strB = '#7E2F8E';
colorB = sscanf(strB(2:end),'%2x%2x%2x',[1 3])/255;

% Self Heating Rate Profile Plot
% Equation 7a plot
subplot (2,4,1)
semilogy(y2a(:, 3), dydt2a(:, 3), 'color', colorA)
hold on
semilogy(y1a(:, 3), dydt1a(:, 3), 'color', colorB)
hold off
xlabel('Temperature (ºC)')
ylabel('log dT/dt (ºC/min)')
xlim([100 220])
ylim([10^-1 10^2.5])
title('Temperature Profile for change in MCMB surface area (Equation A)')
legend('7a a_H', '7a a_L')

% Equation 7b plot
subplot(2,4,2)
semilogy(y2b(:, 3), dydt2b(:, 3))
hold on
semilogy(y1b(:, 3), dydt1b(:, 3))
hold off
xlabel('Temperature (ºC)')
ylabel('log dT/dt (ºC/min)')
xlim([100 220])
ylim([10^-1 10^2.5])
title('Temperature Profile for change in MCMB surface area (Equation B)')
legend('7b a_H', '7b a_L')


% Change in the amount of intercalated Li plot
% Equation 7a plot
subplot(2, 4, 3)
plot(y2a(:, 3), y2a(:, 2), 'color', colorA)
hold on
plot(y1a(:, 3), y1a(:, 2), 'color', colorB)
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of intercalated Li')
xlim([100 220])
ylim([0.4 0.8])
title('Change in Intercalated Li for change in MCMB surface area Equation A)')
legend('7a a_H', '7a a_L')

% Equation 7b plot
subplot(2, 4, 4)
plot(y2b(:, 3), y2b(:, 2))
hold on
plot(y1b(:, 3), y1b(:, 2))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of intercalated Li')
xlim([100 220])
ylim([0.4 0.8])
title('Change in Intercalated Li for change in MCMB surface area (Equation B)')
legend('7b a_H', '7b a_L')

% Change in the amount of metastable SEI plot
% Equation 7a plot
subplot(2, 4, 5)
plot(y2a(:, 3), y2a(:, 1), 'color', colorA)
hold on
plot(y1a(:, 3), y1a(:, 1), 'color', colorB)
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of Metastable SEI')
xlim([100 220])
ylim ([0 0.16])
title('Change in Metastable SEI for change in MCMB surface area (Equation A)')
legend('7a a_H', '7a a_L')

% Equation 7b plot
subplot(2, 4, 6)
plot(y2b(:, 3), y2b(:, 1))
hold on
plot(y1b(:, 3), y1b(:, 1))
hold off
xlabel('Temperature (ºC)')
ylabel('Amount of Metastable SEI')
xlim([100 220])
ylim ([0 0.16])
title('Change in Metastable SEI for change in MCMB surface area (Equation B)')
legend('7b a_H', '7b a_L')


% Growth of the SEI Layer Plot
% Equation 7a plot
subplot(2, 4, 7)
plot(y2a(:, 3), y2a(:, 4), 'color', colorA)
hold on
plot(y1a(:, 3), y1a(:, 4), 'color', colorB)
hold off
xlabel('Temperature (ºC)')
ylabel('Growth of SEI layer')
xlim([100 220])
ylim([0 0.5])
title('Growth of the SEI Layer for change in MCMB surface area (Equation A)')
legend('7a a_H', '7a a_L')

% Equation 7b plot
subplot(2, 4, 8)
plot(y2b(:, 3), y2b(:, 4))
hold on
plot(y1b(:, 3), y1b(:, 4))
hold off
xlabel('Temperature (ºC)')
ylabel('Growth of SEI layer')
xlim([100 220])
ylim([0 0.5])
title('Growth of the SEI Layer for change in MCMB surface area (Equation B)')
legend('7b a_H', '7b a_L')



%% Functions

function dydt = scenario3E7b(t, y, a)
    % Parameter Values
    % Note: All parameter values were determed based on the values found in
    % Richard et all 1999 II and their experiments.
    h1C = 150; % represents h1/C
    h2C = 325; % represents h2/C
    A1 = 1*10^17;
    A2 = 1*10^8;
    E1 = 1.4;
    E2 = 0.8;
    kB = 8.617*10^(-5);
    z0 = 0.15;
    n = 0.5; 
    m = 1.0;
    a0 = 1;

% Defines Differential Equations
    dydt = zeros(4,1);
    
    % Solves for xf 
    dydt(1) = -A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n; 
    
    % Solves for xi
    dydt(2) = -A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*(exp(-y(4)/z0));
    
    % Solves for T
    dydt(3) = h1C*A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n ...
        + h2C*A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*(exp(-y(4)/z0));

    % Solves for z
    dydt(4) = A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(exp(-y(4)/z0));
end

function dydt = scenario3E7a(t, y, a)
    % Parameter Values
    % Note: All parameter values were determed based on the values found in
    % Richard et all 1999 II and their experiments.
    h1C = 150;
    h2C = 325;
    A1 = 1*10^17;
    A2 = 70;
    E1 = 1.4;
    E2 = 0.37;
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
    dydt(2) = -A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*(z0/y(4));
    
    % Solves for T
    dydt(3) = h1C*A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n ...
        + h2C*A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*(z0/y(4));

    % Solves for z
    dydt(4) = A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(z0/y(4));
end






