function [C_bar, tau_bar] = steady_state(params)

    % Calculate and Return Steady State of the social planner problem;

    % Unpack parameters;
    beta       = params.beta;
    alpha      = params.alpha;
    phi        = params.phi;
    gamma      = params.gamma;
    theta_bar  = params.theta_bar;
    A          = params.A;
    xi         = params.xi;

    % Consumption;
    % TODO : Check
    C_bar       = (beta * xi * (1 - phi) * alpha /(1 - beta * phi) + A/theta_bar )^(-1/gamma);           % Eq (12)
    % Taxes;
    % TODO : CHECK
    tau_bar     = (1 + (A * (1-beta*phi)) / alpha * theta_bar * beta * (1-phi))^(-1);                    %TODO: Nuber eq

end

% C;