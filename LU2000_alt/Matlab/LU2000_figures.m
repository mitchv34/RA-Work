clear all
close all
clc

% Parameters
psi = 0.9;
beta = 0.97;
alpha = 0.8;
phi = 0;
gamma = [0.5 1.5];

delta = beta * (phi + alpha * (1 - phi));

% Other parameters?
theta_bar = 1;
A = 1;


params_gamma_low =  create_params(psi, beta, alpha, phi, gamma(1), theta_bar, A);
params_gamma_high = create_params(psi, beta, alpha, phi, gamma(2), theta_bar, A);

% Get Steady State Values
[C_ss_low, tau_ss_low]  = steady_state(params_gamma_low);  
ss_gamma_low = [C_ss_low, tau_ss_low];
[C_ss_high, tau_ss_high]  = steady_state(params_gamma_high);
ss_gamma_high = [C_ss_high, tau_ss_high];

% Dynamics
[C_low, c_low, theta_low, tax_ratio_low] = get_dynamics(params_gamma_low, ss_gamma_low, 1, 100);
[C_high, c_high, theta_high, tax_ratio_high] = get_dynamics(params_gamma_high, ss_gamma_high, 1, 100);

tau_ratio_ss_low = tau_ss_low / ( 1 - tau_ss_low);
tau_ratio_ss_high = tau_ss_high / ( 1 - tau_ss_high);

% Figures
% 
f_low = gen_plots(C_low, c_low, theta_low, tax_ratio_low, params_gamma_low, ss_gamma_low, 21);
set(gca, 'YLim', [0.0 1.4001]);
f_high = gen_plots(C_high, c_high, theta_high, tax_ratio_high, params_gamma_high, ss_gamma_high, 21);
set(gca, 'YLim', [0.0 1.0]);