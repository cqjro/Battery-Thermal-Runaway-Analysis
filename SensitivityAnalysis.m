% CHE210 Project
% Group: 7



% Sensitivity Analysis
% In this section we will be changing values and seeing how a dTdt vs T
% graph changes accordingly using 7a equation
% =========================================================================%
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


h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10;
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
% defining options for the ode solver
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);


%Low Surface Area Sensitivity Analysis
%Changing h1C
tspanA = [0;0.1;1000];
tspanB = [0;0.1;1000];
y0 = [xf0, xi0(1), T0, z0];



h1C = [60, 90, 120, 150, 180]; %Changing h1C values

% solving ode with the different h1C values
[t1b, y1b] = ode15s(@(t, y) scenario3E7b(t, y, aL, h1C(1), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7b(t, y, aL, h1C(2), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7b(t, y, aL, h1C(3), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7b(t, y, aL, h1C(4), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7b(t, y, aL, h1C(5), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);
dydt2b = zeros(length(y2b), 4);
dydt3b = zeros(length(y3b), 4);
dydt4b = zeros(length(y4b), 4);
dydt5b = zeros(length(y5b), 4);


% obtaining self-heating profile using values from ode solver
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7a(t1b(i), y1b(i, :), aL,h1C(1), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7a(t2b(i), y2b(i, :), aL,h1C(2), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7a(t3b(i), y3b(i, :), aL,h1C(3), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7a(t4b(i), y4b(i, :), aL,h1C(4), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7a(t5b(i), y5b(i, :), aL,h1C(5), h2C, A1, A2b, E1, E2b, kB, z0, n, m, a0);
end

%plotting the different h1C self heating rates
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))

hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature Profile with Changing h1C values using equation 7b')
legend('h1C = 60', 'h1C = 90', 'h1C = 120', 'h1C = 150', 'h1C = 180')


%%
%Changing h2C
clear %this resets our y and t values so the ODE function is able to work
h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10; %Were starting with low surface area so we use 0.10
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);
y0 = [xf0, xi0(1), T0, z0];

h2C = [250, 275, 300, 325, 350]; %Changing h2C values

% getting all the values from the ode solvers for the four functions
[t1b, y1b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C(1), A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C(2), A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C(3), A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C(4), A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C(5), A1, A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);

dydt2b = zeros(length(y2b), 4);

dydt3b = zeros(length(y3b), 4);

dydt4b = zeros(length(y4b), 4);

dydt5b = zeros(length(y5b), 4);


%pluggging ode results and getting the exact self-heating rate 
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL,h1C, h2C(1), A1, A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL,h1C, h2C(2), A1, A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL,h1C, h2C(3), A1, A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7b(t4b(i), y4b(i, :), aL,h1C, h2C(4), A1, A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7b(t5b(i), y5b(i, :), aL,h1C, h2C(5), A1, A2b, E1, E2b, kB, z0, n, m, a0);
end


%plotting different h2C temperature profile
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))
hold on
hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature Profile with Changing h2C values using equation 7b')
legend('h2C = 250', 'h2C = 275', 'h2C = 300', 'h2C = 325', 'h2C = 350')
%%

%Changing A1
h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10; %Were starting with low surface area so we use 0.10
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);
y0 = [xf0, xi0(1), T0, z0];

A1 = [1*10^17, 2*10^17, 3*10^17, 4*10^17, 5*10^17]; %Changing h2C values

%obtaining values for the four parameters (xi, xf, z, T) with the ode solver
[t1b, y1b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1(1), A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1(2), A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1(3), A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1(4), A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1(5), A2b, E1, E2b, kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);

dydt2b = zeros(length(y2b), 4);

dydt3b = zeros(length(y3b), 4);

dydt4b = zeros(length(y4b), 4);

dydt5b = zeros(length(y5b), 4);


%getting self-heating rates from the obtained ode results
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL,h1C, h2C, A1(1), A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL,h1C, h2C, A1(2), A2b, E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL,h1C, h2C, A1(3), A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7b(t4b(i), y4b(i, :), aL,h1C, h2C, A1(4), A2b, E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7b(t5b(i), y5b(i, :), aL,h1C, h2C, A1(5), A2b, E1, E2b, kB, z0, n, m, a0);
end

%plotting different A1 temperature profile
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))
hold on
hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature Profile with Changing A1 values using equation 7b')
legend('A1 = 1*10^{17}', 'A1 = 2*10^{17}', 'A1 = 3*10^{17}', 'A1 = 4*10^{17}', 'A1 = 5*10^{17}')

%%
%Changing A2b
clear %this resets our y and t values so the ODE function is able to work
h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10; %Were starting with low surface area so we use 0.10
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);
y0 = [xf0, xi0(1), T0, z0];

A2b = [1*10^8, 2*10^8, 3*10^8, 4*10^8, 5*10^8]; %Changing h2C values

%obtaining the four parameters in the four functions using ode solver
[t1b, y1b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b(1), E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b(2), E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b(3), E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b(4), E1, E2b, kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b(5), E1, E2b, kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);

dydt2b = zeros(length(y2b), 4);

dydt3b = zeros(length(y3b), 4);

dydt4b = zeros(length(y4b), 4);

dydt5b = zeros(length(y5b), 4);


%obtaining self-heating rate using ode values
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL,h1C, h2C, A1, A2b(1), E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL,h1C, h2C, A1, A2b(2), E1, E2b, kB, z0, n, m, a0);
end

for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL,h1C, h2C, A1, A2b(3), E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7b(t4b(i), y4b(i, :), aL,h1C, h2C, A1, A2b(4), E1, E2b, kB, z0, n, m, a0);
end


for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7b(t5b(i), y5b(i, :), aL,h1C, h2C, A1, A2b(5), E1, E2b, kB, z0, n, m, a0);
end


%plotting self heting using different A2
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))
hold on
hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature Profile with Changing A2 values using equation 7b')
legend('A2 = 1*10^8', 'A2 = 2*10^8', 'A2 = 3*10^8', 'A2 = 4*10^8', 'A2 = 5*10^8')

%%
%Changing E1 values
clear %this resets our y and t values so the ODE function is able to work
h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10; %Were starting with low surface area so we use 0.10
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);
y0 = [xf0, xi0(1), T0, z0];

E1 = [1.38, 1.39, 1.4, 1.41, 1.42]; %Changing E1 values

%obtaining values for the four parameters (xi, xf, z, T) with the ode solver
[t1b, y1b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1(1), E2b, kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1(2), E2b, kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1(3), E2b, kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1(4), E2b, kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1(5), E2b, kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);

dydt2b = zeros(length(y2b), 4);

dydt3b = zeros(length(y3b), 4);

dydt4b = zeros(length(y4b), 4);

dydt5b = zeros(length(y5b), 4);


%obtaining self-heating rate using ode values
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL,h1C, h2C, A1, A2b, E1(1), E2b, kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL,h1C, h2C, A1, A2b, E1(2), E2b, kB, z0, n, m, a0);
end

for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL,h1C, h2C, A1, A2b, E1(3), E2b, kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7b(t4b(i), y4b(i, :), aL,h1C, h2C, A1, A2b, E1(4), E2b, kB, z0, n, m, a0);
end


for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7b(t5b(i), y5b(i, :), aL,h1C, h2C, A1, A2b, E1(5), E2b, kB, z0, n, m, a0);
end


%Plotting self heating profile with changing E1
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))
hold on
hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature profile with Changing E1 values using equation 7b')
legend('E1 = 1.38', 'E1 = 1.39', 'E1 = 1.40', 'E1 = 1.41', 'E1 = 1.42')


%Changing E2 values
clear %this resets our y and t values so the ODE function is able to work
h1C = 150;
h2C = 325;
A1 = 1*10^17;
A2b = 1*10^8;
A2a = 70;
E1 = 1.4;
E2b = 0.8;
E2a = 0.37;
kB = 8.617*10^(-5);
z0 = 0.15;
n = 0.5; 
m = 1; 
a0 = 1;
T0 = 100;
xf0 = 0.10; %Were starting with low surface area so we use 0.10
xi0 = [0.75, 0.45, 0.25]; % for 0.0V, 0.089V & 0.127V
z0 = 0.15;
aL = 1.1;
aH = 6.9;
tspan = [0:0.1:1000]; % larger tspan than this results in very long run time
options = odeset('RelTol',1e-6, 'AbsTol', 1e-9, 'NonNegative',[1 2 3 4]);
y0 = [xf0, xi0(1), T0, z0];

E2b = [0.78, 0.8, 0.82, 0.84, 0.86]; %Changing E2 values

%obtaining values for the four parameters (xi, xf, z, T) with the ode solver
[t1b, y1b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1, E2b(1), kB, z0, n, m, a0), tspan, y0, options);
[t2b, y2b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1, E2b(2), kB, z0, n, m, a0), tspan, y0, options);
[t3b, y3b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1, E2b(3), kB, z0, n, m, a0), tspan, y0, options);
[t4b, y4b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1, E2b(4), kB, z0, n, m, a0), tspan, y0, options);
[t5b, y5b] = ode15s(@(t, y) scenario3E7a(t, y, aL, h1C, h2C, A1, A2b, E1, E2b(5), kB, z0, n, m, a0), tspan, y0, options);


dydt1b = zeros(length(y1b), 4);

dydt2b = zeros(length(y2b), 4);

dydt3b = zeros(length(y3b), 4);

dydt4b = zeros(length(y4b), 4);

dydt5b = zeros(length(y5b), 4);


%obtaining self-heating rate using ode values
for i = 1:length(y1b)
   dydt1b(i, :) = scenario3E7b(t1b(i), y1b(i, :), aL,h1C, h2C, A1, A2b, E1, E2b(1), kB, z0, n, m, a0);
end

for i = 1:length(y2b)
   dydt2b(i, :) = scenario3E7b(t2b(i), y2b(i, :), aL,h1C, h2C, A1, A2b, E1, E2b(2), kB, z0, n, m, a0);
end

for i = 1:length(y3b)
   dydt3b(i, :) = scenario3E7b(t3b(i), y3b(i, :), aL,h1C, h2C, A1, A2b, E1, E2b(3), kB, z0, n, m, a0);
end


for i = 1:length(y4b)
   dydt4b(i, :) = scenario3E7b(t4b(i), y4b(i, :), aL,h1C, h2C, A1, A2b, E1, E2b(4), kB, z0, n, m, a0);
end


for i = 1:length(y5b)
   dydt5b(i, :) = scenario3E7b(t5b(i), y5b(i, :), aL,h1C, h2C, A1, A2b, E1, E2b(5), kB, z0, n, m, a0);
end

%plotting self heating rates with changing E2
figure
plot(y1b(:, 3), dydt1b(:, 3))
hold on
plot(y2b(:, 3), dydt2b(:, 3))
hold on
plot(y3b(:, 3), dydt3b(:, 3))
hold on
hold on
plot(y4b(:, 3), dydt4b(:, 3))
hold on
plot(y5b(:, 3), dydt5b(:, 3))
xlim([80 240])
xlabel('Temperature (ºC)')
ylabel('dTdt (ºC/min)')
title('Temperature Profile with Changing E2 values using equation 7b')
legend('E2 = 0.78', 'E2 = 0.8', 'E2 = 0.82', 'E2 = 0.84', 'E2 = 0.86')



function dydt = scenario3E7b(t, y, a, h1C, h2C, A1, A2, E1, E2, kB, z0, n, m, a0)
    

% Defines Differential Equations
    dydt = zeros(4,1);
    
    % Solves for xf 
    dydt(1) = -A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n; 
    
    % Solves for xi
    dydt(2) = -A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*(a/a0)*exp(-y(4)/z0);
    
    % Solves for T
    dydt(3) = h1C*A1*exp(-E1/(kB*(y(3)+273.15)))*abs(y(1))^n ...
        + h2C*A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)*(a/a0)*exp(-y(4)/z0);

    % Solves for z
    dydt(4) = A2*exp(-E2/(kB*(y(3)+273.15)))*y(2)^m*exp(-y(4)/z0);
end

function dydt = scenario3E7a(t, y, a, h1C, h2C, A1, A2, E1, E2, kB, z0, n, m, a0)    
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
