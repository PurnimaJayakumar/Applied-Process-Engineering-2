%% =========================================================
%  ADIABATIC COMPRESSION: H2/N2 MIXTURE
%  Constant gamma assumption: gamma = 1.41
%  W = gamma/(gamma-1) * R * T1 * [(P2/P1)^((gamma-1)/gamma) - 1]
%  Pressure: P1 = 10 bar  -->  P2 = 200 bar
%  Initial Temperature: T1 = 300 K
%  Moles: H2 = 4412 kmol,  N2 = 1470.5 kmol
%% =========================================================

clc; clear; close all;

%% -------------------------------------------------------
%  1. PARAMETERS
% -------------------------------------------------------
R       = 8.314;        % J/(mol.K)
gamma   = 1.41;         % constant gamma (given)
T1      = 300;          % K
P1      = 10;           % bar
P2      = 200;          % bar

n_H2    = 4412e3;       % mol
n_N2    = 1470.5e3;     % mol
n_total = n_H2 + n_N2;  % mol

%% -------------------------------------------------------
%  2. FINAL TEMPERATURE  (isentropic relation)
%     T2/T1 = (P2/P1)^((gamma-1)/gamma)
% -------------------------------------------------------
T2 = T1 * (P2/P1)^((gamma-1)/gamma);

fprintf('Initial Temperature T1 = %.2f K  (%.2f degC)\n', T1, T1-273.15);
fprintf('Final Temperature   T2 = %.2f K  (%.2f degC)\n', T2, T2-273.15);
fprintf('Temperature Rise  DeltaT = %.2f K\n\n', T2-T1);

%% -------------------------------------------------------
%  3. ADIABATIC WORK PER MOLE
%     W/mol = gamma/(gamma-1) * R * T1 * [(P2/P1)^((gamma-1)/gamma) - 1]
%           = Cp * (T2 - T1)          [same thing]
% -------------------------------------------------------
W_per_mol = (gamma/(gamma-1)) * R * T1 * ((P2/P1)^((gamma-1)/gamma) - 1);  % J/mol

% Verify: also = Cp*(T2-T1)
Cp = gamma * R / (gamma - 1);   % J/(mol.K)
W_check = Cp * (T2 - T1);

%% -------------------------------------------------------
%  4. TOTAL WORK
% -------------------------------------------------------
W_total_J  = W_per_mol * n_total;       % J
W_total_kJ = W_total_J / 1e3;           % kJ
W_total_MJ = W_total_kJ / 1e3;          % MJ

%% -------------------------------------------------------
%  5. RESULTS
% -------------------------------------------------------
fprintf('=== PARAMETERS ===\n');
fprintf('gamma              = %.2f\n',   gamma);
fprintf('Cp (mixture)       = %.4f  J/(mol.K)\n', Cp);
fprintf('Cv (mixture)       = %.4f  J/(mol.K)\n', Cp - R);
fprintf('P1 = %g bar,  P2 = %g bar\n\n', P1, P2);

fprintf('=== WORK RESULTS ===\n');
fprintf('W per mole (mixture)     = %.4f  J/mol\n',  W_per_mol);
fprintf('W check via Cp*(T2-T1)   = %.4f  J/mol\n',  W_check);
fprintf('Total moles              = %.1f  kmol\n',   n_total/1e3);
fprintf('Total adiabatic work     = %.4e  kJ\n',     W_total_kJ);
fprintf('Total adiabatic work     = %.4e  MJ\n\n',   W_total_MJ);

fprintf('=== SUMMARY TABLE ===\n');
fprintf('%-32s %15s\n', 'Quantity', 'Value');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-32s %15.2f  K\n',   'T1 (initial)',          T1);
fprintf('%-32s %15.2f  K\n',   'T2 (final)',             T2);
fprintf('%-32s %15.2f  K\n',   'Delta T',                T2-T1);
fprintf('%-32s %15.2f  bar\n', 'P1',                     P1);
fprintf('%-32s %15.2f  bar\n', 'P2',                     P2);
fprintf('%-32s %15.4f\n',      'gamma',                  gamma);
fprintf('%-32s %15.4f  J/mol','W per mole',              W_per_mol);
fprintf('\n');
fprintf('%-32s %15.4e  kJ\n',  'Total W (kJ)',           W_total_kJ);
fprintf('%-32s %15.4e  MJ\n',  'Total W (MJ)',           W_total_MJ);