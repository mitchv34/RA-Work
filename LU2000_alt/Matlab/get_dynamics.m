function [C, c, theta, tau_ratio] = get_dynamics(params, ss, shock, T)


    % Compute the dynamics of the system starting from an x% shock to the
    % produtivyty parameter.
    % Input:
    %   params: structure with the parameters of the model
    %   ss: vector with the steady state values of the model [C, tau]
    %   shock: percentage shock to the productivity parameter
    %   T: number of periods to simulate
    % Output:
    %   C: vector with the consumption dynamics (with tax adjustment)
    %   c: vector with the consumption dynamics (without tax adjustment)
    %   theta: vector with technology shocks
    %   tau_ratio: vector with the tax ratio dynamics



    % Convert from percetage to decimal
    shock = 1 + shock/100; % Eg 1.05 for a 5% shock

    % Unpack parameters
    psi       = params.psi;
    beta      = params.beta;
    alpha     = params.alpha;
    phi       = params.phi;
    gamma     = params.gamma;
    theta_bar = params.theta_bar;
    A         = params.A;
    xi        = params.xi;

    % Unpack the steady state values
    C_ss = ss(1);
    tau_ss = ss(2);

    % Pre allocate 

    theta     = zeros(T + 1, 1); % Productivity
    C         = zeros(T + 1, 1); % Consumption with tax adjustment
    c         = zeros(T + 1, 1); % Consumption without tax adjustment
    X         = zeros(T + 1, 1); % Measure of past average consumption associated with C
    x         = zeros(T + 1, 1); % Measure of past average consumption associated with c
    tau_ratio = zeros(T + 1, 1); % Taxes to after-tax income ratio


    % Set the initial conditions
    C_last = C_ss; % Consumption with tax adjustment is the steady state value
    c_last = C_ss; % Consumption without tax adjustment is the steady state value
    X_last = alpha * C_ss; % Measure of past average consumption associated with C
    x_last = alpha * C_ss; % Measure of past average consumption associated with c
    theta(1) = shock * theta_bar; % Productivity is the steady state value times the shock
    

    for t  = 1:T % Iterate the system  
        
        % Compute current past average consumption
        X(t) = (1 - phi) * alpha * C_last + phi * X_last ;                                  % With tax adjustment
        x(t) = (1 - phi) * alpha * c_last + phi * x_last ;                                  % Without tax adjustment
        
        % Compute current consumption
        C(t) = (beta * xi * (1 - phi) * alpha /(1 - beta * phi) + A/theta(t) )^(-1/gamma);  % With tax adjustment
        c(t) = ((theta(t)/A)*(1-tau_ss))^(1/gamma); % Without tax adjustment
        
        % Compute current tax ratio consistent with the current consumption (C)
        tau_t =  1 - A./theta(t) .* C(t)^(gamma);                                           % Eq 8
        tau_ratio(t) = tau_t / (1 - tau_t );
        % Compute next period productivity
        if t < T
            theta_inv = ( (1 - psi)/theta_bar + psi/theta(t) );
            theta(t+1) = 1/theta_inv;
        end % if
            
        % Update last period values
        X_last = X(t);
        x_last = x(t);
        C_last = C(t);
        c_last = c(t);

    end % End of iteration

end  %get_dynamics