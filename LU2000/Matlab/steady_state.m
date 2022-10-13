function [C_bar, tau_bar] = steady_state(params)

    % Calculate and Return Steady State of the social planner problem;

    % Unpack parameters;
    beta       = params.beta;
    alpha      = params.alpha;
    phi        = params.phi;
    gamma      = params.gamma;
    theta_bar  = params.theta_bar;
    A          = params.A;

    % Consumption;
    C_bar       = (1/(1-alpha))*((theta_bar/A)*(1-(alpha*beta*(1-phi))/(1-beta*phi)))^(1 / gamma);      % Page 360;
    % Taxes;
    tau_bar     = alpha*beta*(1 - phi)/(1-beta*phi);                                                    % Page 360;

end

% C;