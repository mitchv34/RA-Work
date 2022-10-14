function [theta] = simulate_shocks_AR1(params, T, N)

    % Simulate N time series of length T with AR(1) process starting at theta_bar
    % Inputs:
    % - params                : structure with the parameters of the model
    % - T                     : number of periods
    % - N                     : number of simulations
    % Outputs:
    % - theta                 : simulated time series of theta (N x T matrix)

    % Unpack parameters
    theta_bar = params.theta_bar;           % steady state
    sigma     = params.sigma;               % standard deviation of the shock
    % Check if persistence parameter is specified
    if isfield(params, 'rho')
        rho = params.rho;                   % persistence parameter
    else
        rho = 1.0;                          % defaults to white noise around theta_bar
    end

    %{ 
    TODO: Check if shocks are still uniform or if they are now normal, for now 
    assume they are still uniform 
    %}

    % Bounds of Uniform Distribution; 
    l_bound   = -sigma*sqrt(3);
    u_bound   = -l_bound;
    % Generate shocks;
    epsilon   = l_bound + (u_bound - l_bound) * rand(N, T-1);

    % Pre allocate memory
    theta = zeros(N, T);
    theta(:, 1) = theta_bar;
    % Simulate AR(1) process
    for t = 2:T
        theta(:, t) = rho * theta(:, t-1) + epsilon(:, t-1);
    end

end  %simulate_shocks_AR1