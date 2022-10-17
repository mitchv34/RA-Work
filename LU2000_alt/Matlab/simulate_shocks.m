function [theta] = simulate_shocks(params, T, N, which_method)

    % Simulate N time series of length T  process starting at theta_bar
    % Inputs:
    % - params                : structure with the parameters of the model
    % - T                     : number of periods
    % - N                     : number of simulations
    % - which_method          : specify which method to use for the simulation
    % Outputs:
    % - theta                 : simulated time series of theta (N x T matrix)

    % Pass variables to appropriate method
    switch which_method
        case 'paper'
            theta = simulate_shocks_paper(params, T, N);
        case 'AR1'
            theta = simulate_shocks_AR1(params, T, N);
        otherwise
            error('Invalid method for simulation')
    end

end  %simulate_shocks