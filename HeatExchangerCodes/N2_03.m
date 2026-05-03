clc; clear all;
% Cooling of N2 and O2 from 30°C to -170°C
% Flow rates in kmol/day
N2_flow = 1472.06026;
O2_flow = 390.9959761;

% Temperature limits (°C)
T_initial = 30;
T_final = -170;

% Convert to Kelvin
T1 = T_initial + 273.15;  % 303.15 K
T2 = T_final + 273.15;    % 103.15 K

% -----------------------------
% Nitrogen Cp function (100-500 K) using t = T/1000
Cp_N2 = @(T) 28.98641 + 1.853978*(T/1000) + (-9.647459)*(T/1000).^2 + 16.63537*(T/1000).^3 + 0.000117./(T/1000).^2;
Cp_O2 = @(T) 31.32234	 + -20.23531*(T/1000) + 57.86644*(T/1000).^2 + (-36.50624)*(T/1000).^3 + (-0.007374)./(T/1000).^2;


% -----------------------------
% Calculate enthalpy change per kmol using integral
deltaH_N2_per_kmol = integral(Cp_N2, T2, T1);  % kJ/kmol
deltaH_O2_per_kmol = integral(Cp_O2, T2, T1);  % kJ/kmol


% Total enthalpy change (kJ/day)
deltaH_N2 = deltaH_N2_per_kmol * N2_flow;
deltaH_O2 = deltaH_O2_per_kmol * O2_flow;
total_deltaH = deltaH_N2 + deltaH_O2;

% Display results
fprintf('ΔH_N2 = %e kJ/day\n', deltaH_N2);
fprintf('ΔH_O2 = %e kJ/day\n', deltaH_O2);
fprintf('Total ΔH = %e kJ/day\n', total_deltaH);

% Convert to MW
MW = total_deltaH * 1000 / (86400 * 1e6);  % kJ/day -> MW
fprintf('Total cooling duty = %.6f MW\n', MW);