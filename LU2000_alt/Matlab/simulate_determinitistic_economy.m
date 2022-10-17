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
%     phi       = params.phi;
    gamma     = params.gamma;
    theta_Bar = params.theta_bar;
    A         = params.A;
    xi        = params.xi;

    % Unpack steady state values
    C_ss      = ss(1);
    C_ls      = (theta_Bar/A)^(1/gamma);

    % Pre allocate 
%     theta     = theta_Bar * ones(T, 1); % Productivity (Fixed)
    C_sp      = C_ss      * ones(T, 1); % Consumption social planner
    C_ls      = C_ls      * ones(T, 1); % Consumption laissez faire (no taxation)
    % C_ls      = zeros(T, 1);            % Consumption laissez faire (no taxation)
    X_sp      = alpha * C_sp;           % Measure of past average consumption associated with C_ls
    X_ls      = alpha * C_ls;           % Measure of past average consumption associated with C_ls
    % X_ls      = zeros(T, 1);            % Measure of past average consumption associated with C_ls
    n_sp      = C_sp / theta_Bar;       % Labor social planner
    n_ls      = C_ls / theta_Bar;       % Labor laissez faire 
    % n_ls      = zeros(T, 1);            % Labor laissez faire (no taxation)

    % % Set the initial conditions
    % C_ls(1) = C_ss;                     % Initial consumtption is set to be the steady state level
    % X_ls(1) = alpha * C_ss;             % Initial measure of past average consumption is set to be the steady state level
    
    % for t  = 2:T % Iterate the system  
        
    %     % Compute current past average consumption
    %     X_ls(t)     = (1 - phi) * alpha * C_ls(t-1)      + phi * X_ls(t-1)      ; % Laissez faire
        
    %     % Compute current consumption     
    %     C_ls(t)     = X_ls(t) + ((theta(t)/A))^(1/gamma);                         % Laissez faire
        
    %     % Compute current labor supply
    %     n_ls(t)     = C_ls(t)/ theta(t);                                          % Laissez faire
        
    % end % End of iteration

    % % Burn the first T_burn periods
    % C_ls = C_ls(T_burn+1:end);
    % n_ls = n_ls(T_burn+1:end);
    % C_sp = C_sp(T_burn+1:end);
    % n_sp = n_sp(T_burn+1:end);
    % X_ls = X_ls(T_burn+1:end);
    % X_sp = X_sp(T_burn+1:end);
    % theta = theta(T_burn+1:end);

    % sim_results = [C_ls';n_ls';C_sp';n_sp';X_ls';X_sp']; 
    sim_results.C_ls = C_ls';
    sim_results.n_ls = n_ls';
    sim_results.C_sp = C_sp';
    sim_results.n_sp = n_sp';
    sim_results.X_ls = X_ls';
    sim_results.X_sp = X_sp';

end  %simulate_determinitistic_economy