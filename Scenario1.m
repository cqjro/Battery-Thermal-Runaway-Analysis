% CHE210 Project
% Group: 7
%John Larsen and Anastasia Dimov


% Scenario 1
% estimate the self-heating rate profile for different reaction
% orders; the calculation is based on Eq. 1 and 7b, 
% 
% =========================================================================%
%                              Parameters   
%
%  A1 -> frequency factor for the conversion metastable SEI to stable SEI
%        min^-1
%  C -> specific heat of the ARC sample, J/(g K)
%  d0 -> initial SEI thickness, m
%  dt -> time increment, min
%  dT/dt -> self-heating rate, K/min
%  dxf -> change in the amount of metastable SEI, unitless
%  dxi -> change in the amount of intercalated lithium, unitless
%  E1 -> activation energy for the conversion of metastable SEI to stable 
%        SEI, eV
%  h1 -> heat of reaction for the formation of stable SEI from metastable 
%        SEI, J/g
%  kB -> Boltzmann’s constant, eV/K
%  n -> reaction order for the conversion of metastable SEI to stable SEI, 
%       unitless
%=========================================================================%

% Parameter Values
% Note: All parameter values were determed based on the values found in
% Richard et all 1999 II and their experiments.



%setting ode solver and options
tspan = [0 1000];
options = odeset('RelTol',1e-6, 'AbsTol',1e-9, 'NonNegative',[1 2]); 

%defining initial conditions
y0 = [0.1, 80];

n = [0.2, 0.333, 0.5, 1];
A1 = 1 * (10^17);
kB = 8.617 * (10^(-5));
h1C = 150;
E1 = 1.4;


temp_80_n1 = [];
temp_80_n2 = [];
temp_80_n3 = [];
temp_80_n4 = [];
xf_80_n1 = [];
xf_80_n2 = [];
xf_80_n3 = [];
xf_80_n4 = [];

%obtaining obtaining values from ode solver with different n
[t,y] = ode45(@(t,y) scenario1(t, y, n(1), A1, kB, h1C, E1), tspan, y0, options); %ODE Solver
temp_80_n1 = y(:, 2);
xf_80_n1 = y(:,1);
[t,y] = ode45(@(t,y) scenario1(t, y, n(2), A1, kB, h1C, E1), tspan, y0, options); %ODE Solver
temp_80_n2 = y(:, 2);
xf_80_n2 = y(:,1);
[t,y] = ode45(@(t,y) scenario1(t, y, n(3), A1, kB, h1C, E1), tspan, y0, options); %ODE Solver
temp_80_n3 = y(:, 2);
xf_80_n3 = y(:,1);
[t,y] = ode45(@(t,y) scenario1(t, y, n(4), A1, kB, h1C, E1), tspan, y0, options); %ODE Solver
temp_80_n4 = y(:, 2);
xf_80_n4 = y(:,1);


y02 = [0.1, 100];

temp_100_n1 = [];
temp_100_n2 = [];
temp_100_n3 = [];
temp_100_n4 = [];
xf_100_n1 = [];
xf_100_n2 = [];
xf_100_n3 = [];
xf_100_n4 = [];

%obtaining obtaining values from ode solver with different n
[t2, y2] = ode45(@(t2,y2) scenario1(t, y2, n(1), A1, kB, h1C, E1), tspan, y02, options);
temp_100_n1 = y2(:, 2);
xf_100_n1 = y2(:, 1);
[t2, y2] = ode45(@(t2,y2) scenario1(t, y2, n(2), A1, kB, h1C, E1), tspan, y02, options);
temp_100_n2 = y2(:, 2);
xf_100_n2 = y2(:, 1);
[t2, y2] = ode45(@(t2,y2) scenario1(t, y2, n(3), A1, kB, h1C, E1), tspan, y02, options);
temp_100_n3 = y2(:, 2);
xf_100_n3 = y2(:, 1);
[t2, y2] = ode45(@(t2,y2) scenario1(t, y2, n(4), A1, kB, h1C, E1), tspan, y02, options);
temp_100_n4 = y2(:, 2);
xf_100_n4 = y2(:, 1);

y03 = [0.1, 130];
temp_130_n1 = [];
temp_130_n2 = [];
temp_130_n3 = [];
temp_130_n4 = [];
xf_130_n1 = [];
xf_130_n2 = [];
xf_130_n3 = [];
xf_130_n4 = [];

%obtaining obtaining values from ode solver with different n
[t3, y3] = ode45(@(t3,y3) scenario1(t, y3, n(1), A1, kB, h1C, E1), tspan, y03, options);
temp_130_n1 =  y3(:,2);
xf_130_n1 = y3(:, 1);
[t3, y3] = ode45(@(t3,y3) scenario1(t, y3, n(2), A1, kB, h1C, E1), tspan, y03, options);
temp_130_n2 =  y3(:,2);
xf_130_n2 = y3(:, 1);
[t3, y3] = ode45(@(t3,y3) scenario1(t, y3, n(3), A1, kB, h1C, E1), tspan, y03, options);
temp_130_n3 =  y3(:,2);
xf_130_n3 = y3(:, 1);
[t3, y3] = ode45(@(t3,y3) scenario1(t, y3, n(4), A1, kB, h1C, E1), tspan, y03, options);
temp_130_n4 =  y3(:,2);
xf_130_n4 = y3(:, 1);







%we have all of our temps and xfs, now we need to iterate through them and find our
%dy/dt
%our xf and their respective temp arrays are the same size, so we can just
%iterate through them




dydt_80_n1 = [];
dydt_80_n2 = [];
dydt_80_n3 = [];
dydt_80_n4 = [];

dydt_100_n1 = [];
dydt_100_n2 = [];
dydt_100_n3 = [];
dydt_100_n4 = [];

dydt_130_n1 = [];
dydt_130_n2 = [];
dydt_130_n3 = [];
dydt_130_n3 = [];

%solving the self heating rate profile equations with the obtained ode
%values
for i = 1:length(temp_80_n1)
    dydt_80_n1(i) = solver(xf_80_n1(i), temp_80_n1(i), n(1), A1, kB, h1C, E1);
end

for i = 1:length(temp_80_n2)
    dydt_80_n2(i) = solver(xf_80_n2(i), temp_80_n2(i), n(2), A1, kB, h1C, E1);
end

for i = 1:length(temp_80_n3)
    dydt_80_n3(i) = solver(xf_80_n3(i), temp_80_n3(i), n(3), A1, kB, h1C, E1);
end

for i = 1:length(temp_80_n4)
    dydt_80_n4(i) = solver(xf_80_n4(i), temp_80_n4(i), n(4), A1, kB, h1C, E1);
end

for i = 1:length(temp_100_n1)
    dydt_100_n1(i) = solver(xf_100_n1(i), temp_100_n1(i), n(1), A1, kB, h1C, E1);
end

for i = 1:length(temp_100_n2)
    dydt_100_n2(i) = solver(xf_100_n2(i), temp_100_n2(i), n(2), A1, kB, h1C, E1);
end

for i = 1:length(temp_100_n3)
    dydt_100_n3(i) = solver(xf_100_n3(i), temp_100_n3(i), n(3), A1, kB, h1C, E1);
end

for i = 1:length(temp_100_n4)
    dydt_100_n4(i) = solver(xf_100_n4(i), temp_100_n4(i), n(4), A1, kB, h1C, E1);
end

for i = 1:length(temp_130_n1)
    dydt_130_n1(i) = solver(xf_130_n1(i), temp_130_n1(i), n(1), A1, kB, h1C, E1);
end

for i = 1:length(temp_130_n2)
    dydt_130_n2(i) = solver(xf_130_n2(i), temp_130_n2(i), n(2), A1, kB, h1C, E1);
end

for i = 1:length(temp_130_n3)
    dydt_130_n3(i) = solver(xf_130_n3(i), temp_130_n3(i), n(3), A1, kB, h1C, E1);
end

for i = 1:length(temp_130_n4)
    dydt_130_n4(i) = solver(xf_130_n4(i), temp_130_n4(i), n(4), A1, kB, h1C, E1);
end


%plotting the reaction rates temperature profile
semilogy(temp_80_n1, (dydt_80_n1), "b");
hold on
semilogy(temp_80_n2, (dydt_80_n2), "r");
hold on
semilogy(temp_80_n3, (dydt_80_n3), "g");
hold on
semilogy(temp_80_n4, (dydt_80_n4), "m");
hold on
semilogy(temp_100_n1, (dydt_100_n1), "b");
hold on
semilogy(temp_100_n2, (dydt_100_n2), "r");
hold on
semilogy(temp_100_n3, (dydt_100_n3), "g");
hold on
semilogy(temp_100_n4, (dydt_100_n4), "m");
hold on
semilogy(temp_130_n1, (dydt_130_n1), "b");
hold on
semilogy(temp_130_n2, (dydt_130_n2), "r");
hold on
semilogy(temp_130_n3, (dydt_130_n3), "g");
hold on
semilogy(temp_130_n4, (dydt_130_n4), "m");
hold on
ylim([0.001 100])
xlabel('T (Temperature in C)')
ylabel('log dTdt (in C/min)')
title('Rate of Temperature Change in Relation to Temperature')
legend('n = 0.2', 'n = 0.333', 'n = 0.5', 'n = 1')








function dydt = scenario1(t, y, n, A1, kB, h1C, E1)
    %Define Differential Equations
    dydt = zeros(2,1);
    %dydt = [xf, T]

    dydt(1) = -A1*exp(-E1/(kB*(y(2)+273.15)))*y(1)^n; % Solves for xf
    dydt(2) = h1C*A1*exp(-E1/(kB*(y(2)+273.15)))*y(1)^n; % Solves for T
end

function deriv = solver(xf, temp, n, A1, kB, h1C, E1)
     deriv = h1C*A1*exp(-E1/(kB*(temp+273.15)))*xf^n;
end
