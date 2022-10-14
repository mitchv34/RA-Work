function [theta] = simulate_shocks_paper(params, T, N)

    % Simulate N time series of length T  process starting at theta_bar as in LU2000
    % Inputs:
    % - params                : structure with the parameters of the model
    % - T                     : number of periods
    % - N                     : number of simulations
    % Outputs:
    % - theta                 : simulated time series of theta (N x T matrix)

    % Unpack parameters
    psi       = params.psi;                 % persistence of the shock
    theta_bar = params.theta_bar;           % steady state
    sigma     = params.sigma;               % standard deviation of the shock
    
    % Bounds of Uniform Distribution; 
    l_bound   = -sigma*sqrt(3);
    u_bound   = -l_bound;
    % Generate shocks;
    epsilon   = l_bound + (u_bound - l_bound) * rand(N, T-1);

    % Pre allocate memory
    theta = zeros(N, T);
    theta(:, 1) = theta_bar;
    % Simulate stochastic process
    for t  = 2:T % Iterate the system  
        % Compute current productivity
        theta_inv    = ( (1 - psi)./theta_bar + psi./theta(:, t-1) ).*(1 + epsilon(:, t-1));
        theta(:, t)     = 1./theta_inv;
    end

end  %simulate_shocks_paper