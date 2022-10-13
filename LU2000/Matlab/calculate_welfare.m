function [Delta] = calculate_welfare(stochastic_economy, det_economy, params, Nsim)

    % This function calculates the fractional reduction in welfare of the three stochastic economies
    % it uses a deterministic economy as a benchmark
    % Inputs:
    % - stochastic_econmy     : structure with the results of the simulation of the stochastic economy
    % - det_economy           : structure with the results of the simulation of the deterministic economy 
    % - params                : structure with the parameters of the model
    % Outputs:
    % - Delta                 : structure with the fractional reduction in welfare of the three stochastic economies

    % Unpack the parameters
    gamma     = params.gamma;
    beta      = params.beta;
    A         = params.A;

    % Unpack the results of the simulation of the stochastic economy (Consumption, and labor)
%     C_sp_star = stochastic_economy(  Nsim+1:2*Nsim,:);
%     C_sp_bar  = stochastic_economy(2*Nsim+1:3*Nsim,:);
%     C_ls      = stochastic_economy(3*Nsim+1:4*Nsim,:);
%     X_sp_star = stochastic_economy(4*Nsim+1:5*Nsim,:);
%     X_sp_bar  = stochastic_economy(5*Nsim+1:6*Nsim,:);
%     X_ls      = stochastic_economy(6*Nsim+1:7*Nsim,:);
%     n_sp_star = stochastic_economy(7*Nsim+1:8*Nsim,:);
%     n_sp_bar  = stochastic_economy(8*Nsim+1:9*Nsim,:);
%     n_ls      = stochastic_economy(9*Nsim+1:10*Nsim,:);

    C_sp_star = stochastic_economy.C_sp_star;
    C_sp_bar = stochastic_economy.C_sp_bar;
    C_ls = stochastic_economy.C_ls;
    X_sp_star = stochastic_economy.X_sp_star;
    X_sp_bar = stochastic_economy.X_sp_bar;
    X_ls = stochastic_economy.X_ls;
    n_sp_star = stochastic_economy.n_sp_star;
    n_sp_bar = stochastic_economy.n_sp_bar;
    n_ls = stochastic_economy.n_ls;

    T         = size(C_sp_star,2);
    
    % Generate vector of beta^t
    beta_t    = beta.^(0:T-1)';

    % Calculate lifetime utility of the stochastic economies
    U_sp_star = mean( ( ((C_sp_star - X_sp_star).^(1 - gamma) - 1) / (1 - gamma) - A * n_sp_star)*beta_t);
    U_sp_bar  = mean( ( ((C_sp_bar - X_sp_bar).^(1 - gamma) - 1)   / (1 - gamma) - A * n_sp_bar) *beta_t);
    U_ls      = mean( ( ((C_ls - X_ls).^(1 - gamma) - 1)           / (1 - gamma) - A * n_ls)     *beta_t);
    
    % Unpack the results of the simulation of the deterministic economy (Consumption, and labor)
    % C_ls = det_economy(1,:);
    % n_ls = det_economy(2,:);
    % C_sp = det_economy(3,:);
    % n_sp = det_economy(4,:);
    % X_ls = det_economy(5,:);
    % X_sp = det_economy(6,:);
    C_ls = det_economy.C_ls;
    n_ls = det_economy.n_ls;
    C_sp = det_economy.C_sp;
    n_sp = det_economy.n_sp;
    X_ls = det_economy.X_ls;
    X_sp = det_economy.X_sp;

    % Create equations
    eq_1 = @(delta) ((((1-delta).*C_sp - X_sp).^(1-gamma)-1)/(1-gamma) - A * n_sp)*beta_t - U_sp_star;
    eq_2 = @(delta) ((((1-delta).*C_sp - X_sp).^(1-gamma)-1)/(1-gamma) - A * n_sp)*beta_t - U_sp_bar;
    eq_3 = @(delta) ((((1-delta).*C_sp - X_sp).^(1-gamma)-1)/(1-gamma) - A * n_sp)*beta_t - U_ls;
    eq_4 = @(delta) ((((1-delta).*C_ls - X_ls).^(1-gamma)-1)/(1-gamma) - A * n_ls)*beta_t - U_ls;

    % Solve for delta
    delta0 = 0.0;
    delta_sp_star = fzero(eq_1, delta0);
    delta_sp_bar  = fzero(eq_2, delta0);
    delta_sp_ls   = fzero(eq_3, delta0);
    delta_ls_ls   = fzero(eq_4, delta0);

    Delta = [delta_sp_star, delta_sp_bar, delta_sp_ls, delta_ls_ls];
end  %calculate_welfare