function params = create_params(psi, beta, xi, alpha, phi, gamma, theta_bar, A, sigma)

    % Create a structure with the parameters of the model
    % Inputs:
    %   psi:        persistence of productivity
    %   beta:       discount factor
    %   alpha:      capital share
    %   phi:        capital depreciation rate
    %   gamma:      elasticity of substitution
    %   theta_bar:  steady state of the technology
    %   A:          scale parameter of the technology 
    % Optional inputs: 
    %   sigma:      standard deviation of the technology shock
    %   rho:        presistence of the technology shock
    % Output:
    %   params:     structure with the parameters of the model

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
    params.xi       = xi;
    switch nargin
        case 9
            params.sigma    = sigma;
        otherwise
            params.sigma    = 0.00;
    end

end % create_params

% C;


