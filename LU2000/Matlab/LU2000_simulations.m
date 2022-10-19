% Header

clear all
close all
clc


tic


% Set Model parameters
beta            = 0.97;
phi             = 0.0;
gamma_vals      = [1.5 5.0];
alpha_vals      = [0.2 0.5 0.8];


% Parameters that do not affect welfare;
theta_bar       = 1;                    
A               = 1;                   

% Set simulation parameters;
% Reference page 363.
T               = 1000;                 % Number of periods
Tburn           = 100;                  % Number of periods to burn
Nsim            = 1E4;                  % Number of simulations
% Specify the stochastic process for the productivity shock
which_method = 'paper'; % 'paper' or 'AR1'
% Specify stochastic process parameters
sigma_vals      = [0.01 0.04];         % Standard deviation of the shock
psi             = 0.9;                 % For 'paper' method

% sigma_vals      = 0.01;
% psi             = 0.9;
% beta            = 0.97;
% phi             = 0.0;
% gamma_vals      = 1.5;
% alpha_vals      = 0.5;



% Pre allocate space for results
sp_sp_star = zeros(1, numel(gamma_vals) * numel(alpha_vals) * numel(sigma_vals));
sp_sp_bar = zeros(1, numel(gamma_vals) * numel(alpha_vals) * numel(sigma_vals));
sp_ls = zeros(1, numel(gamma_vals) * numel(alpha_vals) * numel(sigma_vals));
ls_det_LS = zeros(1, numel(gamma_vals) * numel(alpha_vals) * numel(sigma_vals));

i = 1;
for sigma = sigma_vals            
    % Set Seed for Identical Simulations Across Parameters;      % Loop over sigma
    % Pre-compute shocks
    rng(90212);              
    disp_params = create_params(psi, 0, 0, 0, 0, theta_bar, 0, sigma); % Create params struct for simulation  
    theta                   = simulate_shocks(disp_params, T + Tburn, Nsim, which_method);
    theta_Bar               = mean(mean(theta)); % Compute mean of theta over simulations and periods
    for gamma = gamma_vals              % Loop over gamma
        for  alpha = alpha_vals         % Loop over alpha
            
            fprintf("Simulating for sigma = " + sigma + " gamma = " + gamma + " alpha = " + alpha + "...")
            
            % Generate Parameters;
            params                  = create_params(psi, beta, alpha, phi, gamma, theta_bar, A, sigma);
            
            % Steady State Calcul ation for Social Planner;
            [C_ss, tau_ss]          = steady_state(params);  
            ss                      = [C_ss, tau_ss];
            
            % Simulate; 
            sim_results             = simulate_economy(params, ss, theta, Tburn);
            % sim_results             = do_simulations(Nsim, params, ss, sigma, T, Tburn);
            
            % Get theta_Bar
%             sim_results_theta       = sim_results(1:Nsim,:);
            % sim_results_theta       = sim_results.theta;

            % Parameters for Deterministic Simulation (needs to include theta_Bar);
            params_det              = params;
            params_det.theta_bar    = theta_Bar;
            
            % Steady State for Deterministic Economy;
            [C_ss_det, tau_ss_det]  = steady_state(params_det);
            ss_det                  = [C_ss_det, tau_ss_det];
            
            % Simulate Deterministic Economy;
            sim_results_det         = simulate_determinitistic_economy(params_det, ss, T);

            % Calculate welfare
            fprintf("Solving for Delta.\n")
            Delta                   = calculate_welfare(sim_results, sim_results_det, params, Nsim);
            
            % Store welfare loss in %
            sp_sp_star(i)   = round( 100 * Delta(1), 3);
            sp_sp_bar(i)    = round( 100 * Delta(2), 3);
            sp_ls(i)        = round( 100 * Delta(3), 3);
            ls_det_LS(i)    = round( 100 * Delta(4), 3);

            i = i + 1;

        end  
    end  
end  

%% Create table

% T = [sp_sp_star,sp_sp_bar,sp_ls,ls_det_LS];
 T = gen_table(sigma_vals, gamma_vals, alpha_vals, sp_sp_star, sp_sp_bar, sp_ls, ls_det_LS);

disp(T)
toc