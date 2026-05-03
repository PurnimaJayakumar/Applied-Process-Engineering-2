clc; clear all;
% Cooling of N2 and O2 from 30°C to -170°C
% Flow rates in kmol/day
N2_flow = 1470.5882;


% Temperature limits (°C)
T_initial = -170;
T_final = 27;

% Convert to Kelvin
T1 = T_initial + 273.15;  % 303.15 K
T2 = T_final + 273.15;    % 103.15 K

% -----------------------------
% Nitrogen Cp function (100-500 K) using t = T/1000
Cp_N2 = @(T) 28.98641 + 1.853978*(T/1000) + (-9.647459)*(T/1000).^2 + 16.63537*(T/1000).^3 + 0.000117./(T/1000).^2;

% -----------------------------
% Calculate enthalpy change per kmol using integral
deltaH_N2_per_kmol = integral(Cp_N2, T1, T2);  % kJ/kmol


% Total enthalpy change (kJ/day)
deltaH_N2 = deltaH_N2_per_kmol * N2_flow;


% Display results
fprintf('ΔH_N2 = %e kJ/day\n', deltaH_N2);


% convert to MW
MW =  deltaH_N2 * 1000 / (86400 * 1e6);  % kJ/day -> MW
fprintf('Total heating duty = %.6f MW\n', MW);