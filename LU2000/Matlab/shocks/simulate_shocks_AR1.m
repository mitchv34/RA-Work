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
    psi       = params.psi;

    % Define AR1 process
    ar1 = arima('Constant',(1-psi) * log(theta_bar),'AR',{psi},'Variance',sigma/10);
    % Simulate shocks
    theta = exp(simulate(ar1,T, 'NumPaths',N)');

end  %simulate_shocks_AR1