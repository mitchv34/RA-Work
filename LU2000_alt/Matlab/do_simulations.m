function [sim_results] = do_simulations(N_sim, params, ss, sigma, T, Tburn)

    %   This function runs the simulations of the economy and returns the results
    %   in the structure sim_results.
    %   Inputs: 
    %   - N_sim: number of simulations
    %   - params: a structure containing the parameters of the model
    %   ss: a vector containing the steady state values of the variables [C_ss, tau_ss]
    %   sigma: a float containing the standard deviation of the shock
    %   T: a float containing the number of periods to simulate
    %   T_burn: a float containing the number of periods to burn (Final sample starts at T_burn+1)
    % Outputs:
    %   sim_results: a structure containing the simulated variables;

    % Preallocate memory
    theta     = zeros(T - Tburn, N_sim);    % Productivity
    C_sp_star = zeros(T - Tburn, N_sim);    % Consumption optimal taxation
    C_sp_bar  = zeros(T - Tburn, N_sim);    % Consumption fixed taxation (steady state level)
    C_ls      = zeros(T - Tburn, N_sim);    % Consumption laissez faire (no taxation)
    X_sp_bar  = zeros(T - Tburn, N_sim);    % Measure of past average consumption associated with C_sp_bar
    X_sp_star = zeros(T - Tburn, N_sim);    % Measure of past average consumption associated with C_sp_star
    X_ls      = zeros(T - Tburn, N_sim);    % Measure of past average consumption associated with C_ls
    n_sp_star = zeros(T - Tburn, N_sim);    % Labor supply under optimal taxation
    n_sp_bar  = zeros(T - Tburn, N_sim);    % Labor supply under fixed taxation
    n_ls      = zeros(T - Tburn, N_sim);    % Labor supply under laissez faire

    for i = 1:N_sim % Simulate economy N_sim times

        sim_results    = simulate_economy(params, ss, sigma, T, Tburn); % Simulate economy
        
        % Store results
        theta(:,i)     = sim_results.theta;
        C_sp_star(:,i) = sim_results.C_sp_star;
        C_sp_bar(:,i)  = sim_results.C_sp_bar;
        C_ls(:,i)      = sim_results.C_ls;
        X_sp_bar(:,i)  = sim_results.X_sp_bar;
        X_sp_star(:,i) = sim_results.X_sp_star;
        X_ls(:,i)      = sim_results.X_ls;
        n_sp_star(:,i) = sim_results.n_sp_star;
        n_sp_bar(:,i)  = sim_results.n_sp_bar;
        n_ls(:,i)      = sim_results.n_ls;

    end  %for

%     sim_results        = [theta';C_sp_star';C_sp_bar';C_ls';X_sp_star';X_sp_bar';X_ls';n_sp_star';n_sp_bar';n_ls'];        
    sim_results.theta     = theta';
    sim_results.C_sp_star = C_sp_star';
    sim_results.C_sp_bar  = C_sp_bar';
    sim_results.C_ls      = C_ls';
    sim_results.X_sp_star = X_sp_star';
    sim_results.X_sp_bar  = X_sp_bar';
    sim_results.X_ls      = X_ls';
    sim_results.n_sp_star = n_sp_star';
    sim_results.n_sp_bar  = n_sp_bar';
    sim_results.n_ls      = n_ls';

    
end  %do_simulations

% C;