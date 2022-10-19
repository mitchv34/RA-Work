function [sim_results] = simulate_determinitistic_economy(params, ss, T)

    % Simulate a deterministic economy
    % Inputs:
    % - params: a structure containing the parameters of the model
    % - ss: a vector containing the steady state values of the variables [C_ss, tau_ss]
    % - T: a float containing the number of periods to simulate
    % - T_burn: a float containing the number of periods to burn (Final sample starts at T_burn+1)
    % Outputs:
    % - sim_results: a structure containing the simulated variables

    % Unpack parameters
    alpha     = params.alpha;
    gamma     = params.gamma;
    theta_Bar = params.theta_bar;
    A         = params.A;

    % Unpack steady state values
    C_ss      = ss(1);
    C_ls      = (theta_Bar/A)^(1/gamma);    

    % Pre allocate 
    C_sp      = C_ss      * ones(T, 1); % Consumption social planner
    C_ls      = C_ls      * ones(T, 1); % Consumption laissez faire (no taxation)
    X_sp      = alpha * C_sp;           % Measure of past average consumption associated with C_ls
    X_ls      = alpha * C_ls;           % Measure of past average consumption associated with C_ls
    n_sp      = C_sp / theta_Bar;       % Labor social planner
    n_ls      = C_ls / theta_Bar;       % Labor laissez faire 

    % Store results
    sim_results.C_ls = C_ls';
    sim_results.n_ls = n_ls';
    sim_results.C_sp = C_sp';
    sim_results.n_sp = n_sp';
    sim_results.X_ls = X_ls';
    sim_results.X_sp = X_sp';

end  %simulate_determinitistic_economy