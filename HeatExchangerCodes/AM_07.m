%% =========================================================
%  ADIABATIC COMPRESSION: H2/N2 MIXTURE (3:1 molar ratio)
%  Cp = A + B*T + C*T^2 + D/T^2  [J/(mol.K)] 
%  Pressure: P1 = 10 bar  -->  P2 = 200 bar
%  Initial Temperature: T1 = 300 K
%  Moles: H2 = 4412 kmol,  N2 = 1470.5 kmol
%  Process: Reversible Adiabatic (Isentropic), dS_mix = 0
% ==========================================================

clc; clear; close all;
R       = 8.314;          % J/(mol.K) universal gas constant
%% -------------------------------------------------------
%  1. SHOMATE COEFFICIENTS  
%     Cp (J/mol.K) = A + B*T + C*T^2 + D/T^2
% -------------------------------------------------------
%            A        B          C           D
coeff_H2 = R*[3.249,  0.422e-3,  0,  0.083e+5];   
coeff_N2 = R*[3.280,  0.593e-3,  0,   0.040e+5];
coeff_NH3 = R*[3.578,  3.020e-3,  0,   -0.186e+5];
% Unpack
A_H2 = coeff_H2(1); B_H2 = coeff_H2(2); C_H2 = coeff_H2(3); D_H2 = coeff_H2(4);
A_N2 = coeff_N2(1); B_N2 = coeff_N2(2); C_N2 = coeff_N2(3); D_N2 = coeff_N2(4);
A_NH3= coeff_NH3(1); B_NH3 = coeff_NH3(2); C_NH3 = coeff_NH3(3); D_NH3 = coeff_NH3(4);

%% -------------------------------------------------------
%  2. SYSTEM PARAMETERS
% -------------------------------------------------------


n_H2    = 17647.0584e3;      
n_N2    = 5882.3528e3;      
n_NH3   = 2.941e3; 
n_total = n_H2 + n_N2 + n_NH3;   % total moles

y_H2    = n_H2 / n_total; % mole fraction H2
y_N2    = n_N2 / n_total; % mole fraction N2
y_NH3   = n_NH3 / n_total;
T1 = -86.3+273;                 % K  initial temperature
P1 = 80;                  % bar initial pressure
P2 = 210;                 % bar final pressure

fprintf('=== SYSTEM SETUP ===\n');
fprintf('n_H2    = %E kmol/day\n', n_H2/1e3);
fprintf('n_N2    = %E kmol/day\n', n_N2/1e3);
fprintf('n_NH3    = %E kmol/day\n', n_NH3/1e3);
fprintf('T1 = %g K,  P1 = %g bar,  P2 = %g bar\n\n', T1, P1, P2);

%% -------------------------------------------------------
%  3. Cp FUNCTIONS  [J/(mol.K)]
% -------------------------------------------------------
Cp_H2 = @(T) A_H2 + B_H2*T + C_H2*T.^2 + D_H2./T.^2;
Cp_N2 = @(T) A_N2 + B_N2*T + C_N2*T.^2 + D_N2./T.^2;
Cp_NH3 = @(T) A_NH3 + B_NH3*T + C_NH3*T.^2 + D_NH3./T.^2;
% Mixture Cp per mole of mixture
Cp_mix = @(T) y_H2 * Cp_H2(T) + y_N2 * Cp_N2(T)+y_NH3 * Cp_NH3(T);

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
dH_H2_per_mol = integral(@(T) Cp_H2(T), T1, T2);  % J/mol H2
dH_N2_per_mol = integral(@(T) Cp_N2(T), T1, T2);  % J/mol N2
dH_NH3_per_mol = integral(@(T) Cp_NH3(T), T1, T2);  % J/mol NH3

W_H2_kJ = dH_H2_per_mol * n_H2 / 1e3;
W_N2_kJ = dH_N2_per_mol * n_N2 / 1e3;
W_NH3_kJ = dH_NH3_per_mol * n_NH3 / 1e3;


fprintf('=== ENTHALPY / WORK RESULTS ===\n');
fprintf('DeltaH mixture    = %.4f  J/mol of mix\n', dH_per_mol);
fprintf('DeltaH H2         = %.4f  J/mol H2\n',     dH_H2_per_mol);
fprintf('DeltaH N2         = %.4f  J/mol N2\n\n',   dH_N2_per_mol);

fprintf('Adiabatic work (H2 contribution) = %.4e  kJ\n', W_H2_kJ);
fprintf('Adiabatic work (N2 contribution) = %.4e  kJ\n', W_N2_kJ);
fprintf('Adiabatic work (NH3 contribution) = %.4e  kJ\n', W_NH3_kJ);
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