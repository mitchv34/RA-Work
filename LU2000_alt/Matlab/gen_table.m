function [Table] = gen_table(sigma_vals, gamma_vals, alpha_vals, sp_sp_star, sp_sp_bar, sp_ls, ls_det_LS)

    % Generate a table with the results of the simulations for the different values of sigma, gamma and alpha
    % Inputs:
    %   sigma_vals: (vector or number) with the values of sigma
    %   gamma_vals: (vector or number) with the values of gamma
    %   alpha_vals: (vector or number) with the values of alpha
    %   sigma_vals, gamma_vals, alpha_vals, sp_sp_star, sp_sp_bar, sp_ls, ls_det_LS: results of the simulations
    % Outputs:
    %   Table: table with the results of the simulations for the different values of sigma, gamma and alpha
    
    
    % Pre allocate columns with the values of the parameters
    sigma_row = zeros(1, length(sigma_vals) * length(gamma_vals) * length(alpha_vals));
    gamma_row = zeros(1, length(sigma_vals) * length(gamma_vals) * length(alpha_vals));
    alpha_row = zeros(1, length(sigma_vals) * length(gamma_vals) * length(alpha_vals));
    i = 1;                                             % iteration counter
    for sigma_i = 1:numel(sigma_vals)                  % Loop over sigma values
        for gamma_i = 1:numel(gamma_vals)              % Loop over gamma values
            for  alpha_i = 1:numel(alpha_vals)         % Loop over alpha values
                sigma_row(i) = sigma_vals(sigma_i);   % Fill in the values of sigma
                gamma_row(i) = gamma_vals(gamma_i);   % Fill in the values of gamma
                alpha_row(i) = alpha_vals(alpha_i);   % Fill in the values of alpha
                i = i + 1;                             % Update the iteration counter
            end
        end
    end

    col_names = ["\sigma", "\gamma", "\alpha", "SP -> SP_star", "SP -> SP_bar", "SP -> LS", "LS_det -> LS"];

    Table = table(sigma_row', gamma_row', alpha_row', sp_sp_star', sp_sp_bar', sp_ls', ls_det_LS', 'VariableNames', col_names);

end  %gen_table