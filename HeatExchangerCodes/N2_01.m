%% =========================================================
%  ADIABATIC COMPRESSION: H2/N2 MIXTURE (3:1 molar ratio)
%  Cp = A + B*T + C*T^2 + D/T^2  [J/(mol.K)] 
%  Pressure: P1 = 1 bar  -->  P2 = 6 bar
%  Initial Temperature: T1 = 300 K
% 
%  Process: Reversible Adiabatic (Isentropic), dS_mix = 0
% ==========================================================

clc; clear; close all;
R       = 8.314;          % J/(mol.K) universal gas constant
%% -------------------------------------------------------
%  1.  COEFFICIENTS  
%     Cp (J/mol.K) = A + B*T + C*T^2 + D/T^2
% -------------------------------------------------------
%            A        B          C           D
coeff_O2 = R*[3.639,  0.506e-3,  0,  -0.227e+5];   
coeff_N2 = R*[3.280,  0.593e-3,  0,   0.040e+5];
coeff_CO2 = R*[5.457,  1.045e-3,  0,   -1.157e+5];
coeff_H2O = R*[3.470,  1.450e-3,  0,   0.121e+5];
% Unpack
A_O2 = coeff_O2(1); B_O2 = coeff_O2(2); C_O2 = coeff_O2(3); D_O2 = coeff_O2(4);
A_N2 = coeff_N2(1); B_N2 = coeff_N2(2); C_N2 = coeff_N2(3); D_N2 = coeff_N2(4);
A_CO2 = coeff_CO2(1); B_CO2 = coeff_CO2(2); C_CO2 = coeff_CO2(3); D_CO2 = coeff_CO2(4);
A_H2O = coeff_H2O(1); B_H2O = coeff_H2O(2); C_H2O = coeff_H2O(3); D_H2O = coeff_H2O(4);

%% -------------------------------------------------------
%  2. SYSTEM PARAMETERS
% -------------------------------------------------------


n_N2    = 1472.06026e3;      
n_O2    = 390.9959761e3;     
n_CO2   = 0.746549694e3;
n_H2O   = 35.25373555e3; 

n_total = n_O2 + n_N2 + n_CO2 + n_H2O; % total moles

y_O2    = n_O2 / n_total; % mole fraction H2
y_N2    = n_N2 / n_total; % mole fraction N2
y_CO2    = n_CO2 / n_total;
y_H2O    = n_H2O / n_total;

T1 = 300;                 % K  initial temperature
P1 = 1;                  % bar initial pressure
P2 = 6;                 % bar final pressure



%% -------------------------------------------------------
%  3. Cp FUNCTIONS  [J/(mol.K)]
% -------------------------------------------------------
Cp_O2 = @(T) A_O2 + B_O2*T + C_O2*T.^2 + D_O2./T.^2;
Cp_N2 = @(T) A_N2 + B_N2*T + C_N2*T.^2 + D_N2./T.^2;
Cp_CO2 = @(T) A_CO2 + B_CO2*T + C_CO2*T.^2 + D_CO2./T.^2;
Cp_H2O = @(T) A_H2O + B_H2O*T + C_H2O*T.^2 + D_H2O./T.^2;

% Mixture Cp per mole of mixture
Cp_mix = @(T) y_O2 * Cp_O2(T) + ...
              y_N2 * Cp_N2(T) + ...
              y_CO2 * Cp_CO2(T) + ...
              y_H2O * Cp_H2O(T);
%% -------------------------------------------------------
%  4. ISENTROPIC CONDITION
%     For ideal gas mixture, isentropic means dS = 0:
%
%     integral(T1 -> T2) [Cp_mix(T)/T] dT  =  R * ln(P2/P1)
%
%     Solve for T2 numerically using fzero
% -------------------------------------------------------

% Left-hand side: entropy change due to temperature
% Right-hand side: entropy change due to pressure (must equal LHS)
RHS = R * log(P2 / P1);   % J/(mol.K)

fprintf('=== ISENTROPIC CONDITION ===\n');
fprintf('R * ln(P2/P1) = %.4f J/(mol.K)\n\n', RHS);

% Define residual function: integral(Cp/T dT) from T1 to T2 - RHS = 0
residual = @(T2) integral(@(T) Cp_mix(T)./T, T1, T2) - RHS;

% Initial guess for T2 (use ideal gamma ~ 1.4 estimate)
T2_guess = T1 * (P2/P1)^((1.4-1)/1.4);

% Solve
options = optimset('TolX', 1e-6, 'Display', 'off');
T2 = fzero(residual, T2_guess, options);

fprintf('=== TEMPERATURE CHANGE ===\n');
fprintf('T1 (initial)   = %.2f K  (%.2f degC)\n', T1, T1-273.15);
fprintf('T2 (final)     = %.2f K  (%.2f degC)\n', T2, T2-273.15);
fprintf('Delta T        = %.2f K\n\n', T2 - T1);

%% -------------------------------------------------------
%  5. ENTHALPY CHANGE = ADIABATIC WORK (reversible, W_s = DeltaH)
%     DeltaH per mole of mixture = integral(T1->T2) Cp_mix(T) dT
% -------------------------------------------------------

dH_per_mol = integral(@(T) Cp_mix(T), T1, T2);   % J/mol of mixture

% Total work  [convert J -> kJ]
W_total_kJ = dH_per_mol * n_total / 1e3;          % kJ

% Per-component breakdown (useful for verification)
dH_O2_per_mol = integral(@(T) Cp_O2(T), T1, T2);  % J/mol O2
dH_N2_per_mol = integral(@(T) Cp_N2(T), T1, T2);  % J/mol N2
dH_CO2_per_mol = integral(@(T) Cp_CO2(T), T1, T2);
dH_H2O_per_mol = integral(@(T) Cp_H2O(T), T1, T2);

W_O2_kJ = dH_O2_per_mol * n_O2 / 1e3;
W_N2_kJ = dH_N2_per_mol * n_N2 / 1e3;
W_CO2_kJ = dH_CO2_per_mol * n_CO2 / 1e3;
W_H2O_kJ = dH_H2O_per_mol * n_H2O / 1e3;

fprintf('=== ENTHALPY / WORK RESULTS ===\n');
fprintf('DeltaH mixture    = %.4f  J/mol of mix\n', dH_per_mol);
fprintf('DeltaH O2         = %.4f  J/mol O2\n',     dH_O2_per_mol);
fprintf('DeltaH N2         = %.4f  J/mol N2\n\n',   dH_N2_per_mol);
fprintf('DeltaH CO2         = %.4f  J/mol CO2\n',     dH_CO2_per_mol);
fprintf('DeltaH H2O         = %.4f  J/mol H2O\n\n',   dH_H2O_per_mol);

fprintf('Adiabatic work (O2 contribution) = %.4e  kJ\n', W_O2_kJ);
fprintf('Adiabatic work (N2 contribution) = %.4e  kJ\n', W_N2_kJ);
fprintf('Adiabatic work (CO2 contribution) = %.4e  kJ\n', W_CO2_kJ);
fprintf('Adiabatic work (H2O contribution) = %.4e  kJ\n', W_H2O_kJ);

fprintf('TOTAL adiabatic work             = %.4e  kJ\n\n', W_total_kJ);


%  7. GAMMA (EFFECTIVE) AT INITIAL AND FINAL STATE
%     gamma = Cp / Cv,  Cv = Cp - R  (ideal gas)
% -------------------------------------------------------
Cp_mix_T1 = Cp_mix(T1);
Cp_mix_T2 = Cp_mix(T2);
Cv_mix_T1 = Cp_mix_T1 - R;
Cv_mix_T2 = Cp_mix_T2 - R;
gamma_T1   = Cp_mix_T1 / Cv_mix_T1;
gamma_T2   = Cp_mix_T2 / Cv_mix_T2;

fprintf('=== EFFECTIVE GAMMA ===\n');
fprintf('Cp_mix at T1 = %.4f J/(mol.K),  gamma at T1 = %.4f\n', Cp_mix_T1, gamma_T1);
fprintf('Cp_mix at T2 = %.4f J/(mol.K),  gamma at T2 = %.4f\n', Cp_mix_T2, gamma_T2);
fprintf('(Note: gamma varies with T, which is why constant-gamma formula is NOT used)\n\n');

%% -------------------------------------------------------


fprintf('=== SUMMARY TABLE ===\n');
fprintf('%-30s %15s\n', 'Quantity', 'Value');
fprintf('%s\n', repmat('-',1,47));
fprintf('%-30s %15.2f  K\n',    'Initial Temperature T1',    T1);
fprintf('%-30s %15.2f  K\n',    'Final Temperature T2',      T2);
fprintf('%-30s %15.2f  K\n',    'Temperature Rise DeltaT',   T2-T1);
fprintf('%-30s %15.2f  bar\n',  'Initial Pressure P1',       P1);
fprintf('%-30s %15.2f  bar\n',  'Final Pressure P2',         P2);
fprintf('%-30s %15.4e  kJ\n',   'Total Adiabatic Work',      W_total_kJ);
fprintf('%-30s %15.4f  J/mol\n','DeltaH per mol mixture',    dH_per_mol);
fprintf('%-30s %15.4f\n',       'gamma at T1',               gamma_T1);
fprintf('%-30s %15.4f\n',       'gamma at T2',               gamma_T2);
fprintf('n_total = %.4e mol\n', n_total);