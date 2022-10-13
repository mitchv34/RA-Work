function params = create_params(psi, beta, alpha, phi, gamma, theta_bar, A)

    % Create a structure with the parameters of the model
    
    % Calculate delta
    delta = beta * (phi + alpha * (1 - phi));

    % Create the structure
    params.psi      = psi;
    params.beta     = beta;
    params.alpha    = alpha;
    params.phi      = phi;
    params.gamma    = gamma;
    params.theta_bar= theta_bar;
    params.A        = A;
    params.delta    = delta;

end % create_params

% C;


