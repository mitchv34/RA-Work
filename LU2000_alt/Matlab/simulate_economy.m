function [sim_results] = simulate_economy(params, ss, theta, T_burn)

    % Simulate economy
    % Inputs:
    % - params: a structure containing the parameters of the model
    % - ss: a vector containing the steady state values of the variables [C_ss, tau_ss]
    % - theta: a a vector containing the simulation of the shock
    % - T_burn: a float containing the number of periods to burn (Final sample starts at T_burn+1)
    % Outputs:
    % - sim_results: a structure containing the simulated variables

    % Unpack parameters
    beta      = params.beta;
    alpha     = params.alpha;
    phi       = params.phi;
    gamma     = params.gamma;
    theta_bar = params.theta_bar;
    A         = params.A;
    xi        = params.xi;

    % Unpack Steady state values
    C_ss      = ss(1);
    tau_ss    = ss(2); 

    % Infer T and N from theta
    [N, T]    = size(theta);

    % Pre-allocate; 
    C_sp_star = zeros(N, T);        % Consumption optimal taxation;
    C_sp_bar  = zeros(N, T);        % Consumption fixed taxation (steady state level);
    C_ls      = zeros(N, T);        % Consumption laissez faire (no taxation);
    X_sp_bar  = zeros(N, T);        % Measure of past average consumption associated with C_sp_bar;
    X_sp_star = zeros(N, T);        % Measure of past average consumption associated with C_sp_star;
    X_ls      = zeros(N, T);        % Measure of past average consumption associated with C_ls;
    n_sp_star = zeros(N, T);        % Labor supply under optimal taxation;
    n_sp_bar  = zeros(N, T);        % Labor supply under fixed taxation;
    n_ls      = zeros(N, T);        % Labor supply under laissez faire;

    % Set the initial conditions
    C_sp_star(:, 1) = C_ss;            % Initial consumtption is steady state level;
    C_sp_bar(:, 1)  = C_ss;           
    C_ls(:, 1)      = 1/(1-alpha)*(theta_bar/A)^(1/gamma);   % Different steady state;         
    X_sp_star(:, 1) = alpha * C_ss;    % Initial measure of past average consumption is steady state level;
    X_sp_bar(:, 1)  = alpha * C_ss;    
    X_ls(:, 1)      = alpha * C_ls(1);    

    for t  = 2:T % Iterate the system  
        
        % Compute current past average consumption
        X_sp_star(:, t) = (1 - phi) * alpha * C_sp_star(:, t-1) + phi * X_sp_star(:, t-1);                     % Optimal;
        X_sp_bar(:, t)  = (1 - phi) * alpha * C_sp_bar(:, t-1)  + phi * X_sp_bar(:, t-1);                      % Fixed tax;
        X_ls(:, t)      = (1 - phi) * alpha * C_ls(:, t-1)      + phi * X_ls(:, t-1);                          % Laissez faire;
        
        % Compute current consumption
        C_sp_star(:, t) = ( ( beta * xi* alpha * (1 - phi) )/( 1- beta * phi) + (A./theta(:, t)) ).^(-1/gamma);  % Optimal; 
        C_sp_bar(:, t)  = ( ( theta(:, t) ./ A)*(1-tau_ss)).^(1./gamma);                                        % Fixed tax;        
        C_ls(:, t)      = ( ( theta(:, t) ./ A)).^(1./gamma);                                                   % Laissez faire;
        
        % Compute current labor supply
        n_sp_star(:, t) = C_sp_star(:, t) ./ theta(:, t);                                                       % Optimal;
        n_sp_bar(:, t)  = C_sp_bar(:, t)  ./ theta(:, t);                                                       % Fixed tax;
        n_ls(:, t)      = C_ls(:, t)      ./ theta(:, t);                                                       % Laissez faire;

    end % End of iteration

    % Burn the first T_burn periods
    C_sp_star = C_sp_star(:, T_burn+1:end);
    C_sp_bar  = C_sp_bar(:, T_burn+1:end);
    C_ls      = C_ls(:, T_burn+1:end);
    X_sp_star = X_sp_star(:, T_burn+1:end);
    X_sp_bar  = X_sp_bar(:, T_burn+1:end);
    X_ls      = X_ls(:, T_burn+1:end);
    n_sp_star = n_sp_star(:, T_burn+1:end);
    n_sp_bar  = n_sp_bar(:, T_burn+1:end);
    n_ls      = n_ls(:, T_burn+1:end);

    % Store results
    sim_results.C_sp_star = C_sp_star';
    sim_results.C_sp_bar  = C_sp_bar';
    sim_results.C_ls      = C_ls';
    sim_results.X_sp_star = X_sp_star';
    sim_results.X_sp_bar  = X_sp_bar';
    sim_results.X_ls      = X_ls';
    sim_results.n_sp_star = n_sp_star';
    sim_results.n_sp_bar  = n_sp_bar';
    sim_results.n_ls      = n_ls';
    
end  %simulate_economy