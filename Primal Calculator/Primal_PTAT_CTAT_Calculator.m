%% PTAT - CTAT Design
% Computes PTAT and CTAT characteristics for a BJT junction.
% Modes:
%   1) Single-temperature point
%   2) Temperature interval sweep
% Units:
%   T in °C, V in V, I in A, energies in J.

clear; clc;

% -- Constants --
k     = 1.380649e-23;    % Boltzmann constant [J/K]
q     = 1.602176634e-19; % Elementary charge [C]
A     = 7.921e-5;        % Junction area [cm^2]
T_0   = 300;             % Reference T [K]

% Recombination parameters
sigma_p = 5.2e-17;       % Hole capture [cm^2]
sigma_n = 5.2e-17;       % Electron capture [cm^2]
v_thp   = 1.65e5;        % Hole thermal velocity [cm/s]
v_thn   = 2.3e5;         % Electron thermal velocity [cm/s]
N_t     = 1e13;          % Trap density [cm^-3]

% Doping
N_A = 1.48e15;           % Acceptor [cm^-3]
N_D = 1.75e16;           % Donor    [cm^-3]

% Mobility exponents
m_p   = -2.4;            m_n   = -2.2;
mu_p0 = 480;             mu_n0 = 1350;        % [cm^2/V·s]

% Bandgap parameters (Si)
E_g0_J = 1.17 * q;       % Bandgap at 0K [J]
alpha  = 4.73e-4 * q;    % Varshni alpha [J/K]
beta   = 636;            % Varshni beta [K]

% Intrinsic prefactor and Auger coeffs
B     = 1.08e31;         % [cm^-6·K^-3]
C_p   = 9.9e-32;         % Auger hole [cm^6/s]
C_n   = 2.8e-31;         % Auger electron [cm^6/s]

mode = input('Point (1) or interval (2) for temperature? ');
assert(ismember(mode,[1,2]), 'Mode must be 1 or 2');

if mode==1
  %% Single-point mode
  T   = input('Enter temperature (°C): ');
  T_k = T + 273;
  Vt  = k * T_k / q;
  fprintf('At T=%.2f°C, V_T=%.4f V\n', T, Vt);

  [tau_p, tau_n] = compute_lifetimes(T_k, sigma_p,v_thp, sigma_n,v_thn, N_t, C_p,C_n, N_A, N_D, B, E_g0_J, k);
  [D_p, D_n]    = compute_transport(T_k, mu_p0,mu_n0, m_p,m_n, k,q);
  [Mp, Mn, M]   = compute_diff_life_factors(D_p,D_n, tau_p,tau_n, N_D, N_A);
  m_diff        = compute_mdiff(Mp, Mn, m_p, m_n);

  I_s = compute_Is(T_k, E_g0_J, alpha, beta, k, B, A, q, M);
  choice = input('Supply V_d (1) or I_c (2)? ');
  assert(ismember(choice,[1,2]),'Choice must be 1 or 2');
  if choice==1
    V_d = input('Enter V_d [V]: ');
  else
    I_c = input('Enter I_c [A]: ');
    V_d = Vt * log(I_c / I_s);
  end

  % CTAT slope
  TC = (V_d - ((4 + m_diff)*Vt) - (E_g0_J/q)) / T_k;
  fprintf('TC = %.4f V/K\n', TC);

  % Solid‐PTAT: ask number of diodes
  n = input('Number of parallel diodes (n)? ');
  I_d = I_s * exp(V_d / Vt);
  I_tot = n * I_d;
  V_d1  = Vt * log(I_tot / I_s);
  PTAT = V_d1 - V_d;      %  ensure positive ΔV
  fprintf('PTAT ΔV = %.4f V\n', PTAT);

  %--- Plot 1: Temperature vs TC ---
  figure;
  plot(T, TC, 'o','LineWidth',1.5);
  xlabel('Temperature (°C)'); ylabel('TC [V/K]');
  title('T vs TC (single point)');

  %--- Plot 2: Temperature vs V_d ---
  figure;
  plot(T, V_d, 'o','LineWidth',1.5);
  xlabel('Temperature (°C)'); ylabel('V_d [V]');
  title('T vs V_d (single point)');

  %--- Plot 3: Temperature vs ΔV (PTAT) ---
  figure;
  plot(T, PTAT, 'o','LineWidth',1.5);
  xlabel('Temperature (°C)'); ylabel('ΔV (V)');
  title('T vs PTAT (V_d1 - V_d) (single point)');

else
  %% Interval mode
  T_start = input('Start T [°C]: ');
  T_end   = input('End   T [°C]: ');
  T_step  = input('Step  [°C]: ');
  assert(T_step>0,'Step must be positive');
  T       = T_start:T_step:T_end;
  T_k     = T + 273;
  Vt      = k .* T_k ./ q;

  [tau_p, tau_n] = compute_lifetimes(T_k, sigma_p,v_thp, sigma_n,v_thn, N_t, C_p,C_n, N_A, N_D, B, E_g0_J, k);
  [D_p, D_n]    = compute_transport(T_k, mu_p0,mu_n0, m_p,m_n, k,q);
  [Mp, Mn, M]   = compute_diff_life_factors(D_p,D_n, tau_p,tau_n, N_D, N_A);
  m_diff        = compute_mdiff(Mp, Mn, m_p, m_n);

  I_s_vec = compute_Is(T_k, E_g0_J, alpha, beta, k, B, A, q, M);
  choice   = input('Supply V_d vector (1) or I_c vector (2)? ');
  assert(ismember(choice,[1,2]),'Choice must be 1 or 2');
  if choice==1
    V_d = input('Enter V_d vector [V]: ');
    assert(numel(V_d)==numel(T),'length mismatch');
    I_c = I_s_vec .* exp(V_d ./ Vt);
  else
    I_c = input('Enter I_c vector [A]: ');
    assert(numel(I_c)==numel(T),'length mismatch');
    V_d = Vt .* log(I_c ./ I_s_vec);
  end

  % CTAT vs T
  TC = (V_d - ((4 + m_diff).*Vt) - (E_g0_J/q)) ./ T_k;
  fprintf('Average TC = %.4f V/K\n', mean(TC));

  % Solid‐PTAT: ask diodes
  n = input('Number of parallel diodes (n)? ');
  I_d = I_s_vec .* exp(V_d ./ Vt);
  I_tot = n .* I_d;
  V_d1  = Vt .* log(I_tot ./ I_s_vec);
  PTAT = V_d1 - V_d;      %  ensure positive ΔV
  fprintf('PTAT ΔV = %.4f V\n', PTAT);

  %--- Plot 1: T vs TC ---
  figure;
  plot(T, TC, '-o','LineWidth',1.5);
  xlabel('Temperature (°C)'); ylabel('TC [V/K]');
  title('T vs TC (interval)');

  %--- Plot 2: T vs V_d & PTAT ---
  figure;
  plot(T, V_d, '-s','LineWidth',1.5); hold on;
  plot(T, PTAT, '-d','LineWidth',1.5);
  xlabel('Temperature (°C)'); ylabel('Voltage [V]');
  title('T vs V_d & V_d1 - V_d (interval)');
  legend('V_d','ΔV (PTAT)','Location','best');
end

%% Local functions
function [tau_p, tau_n] = compute_lifetimes(Tk, sp,vp, sn,vn, Nt, Cp,Cn, NA,ND, B, Eg0, k)
  tau_SRHp = 1/(sp*vp*Nt);
  tau_SRHn = 1/(sn*vn*Nt);
  tau_Aupr = 1/(Cp*NA^2);
  tau_Aun  = 1/(Cn*ND^2);
  ni       = sqrt(B*Tk.^3 .* exp(-Eg0./(k*Tk)));
  tau_rad  = 1./(4.7e-15*2.*ni);
  tau_p    = 1./(1/tau_SRHp + 1/tau_Aupr + 1./tau_rad);
  tau_n    = 1./(1/tau_SRHn + 1/tau_Aun + 1./tau_rad);
end

function [Dp, Dn] = compute_transport(Tk, mup0,mun0, mp, mn, k, q)
  mup = mup0 .* (Tk/300).^mp;
  mun = mun0 .* (Tk/300).^mn;
  Dp  = mup .* (k.*Tk./q);
  Dn  = mun .* (k.*Tk./q);
end

function [Mp, Mn, M] = compute_diff_life_factors(Dp, Dn, taup, taun, ND, NA)
  Mp = (1/ND) * sqrt(Dp./taup);
  Mn = (1/NA) * sqrt(Dn./taun);
  M  = Mp + Mn;
end

function m = compute_mdiff(Mp, Mn, mp, mn)
  Msum = Mp + Mn;
  m    = ((mp+1)/2 .* Mp + (mn+1)/2 .* Mn) ./ Msum;
end

function Is = compute_Is(Tk, Eg0, alpha, beta, k, B, A, q, M)
  Eg   = Eg0 - alpha.*Tk.^2 ./ (Tk + beta);
  ni2  = B .* Tk.^3 .* exp(-Eg./(k.*Tk));
  Is   = q * A .* ni2 .* M;
end
