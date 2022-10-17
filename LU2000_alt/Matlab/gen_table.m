function [Table] = gen_table(sigma_vals, gamma_vals, alpha_vals, sp_sp_star, sp_sp_bar, sp_ls, ls_det_LS)


    % Generate a table with the results of the simulations


    sigma_row =[sigma_vals(1) sigma_vals(1) sigma_vals(1) sigma_vals(1) sigma_vals(1) sigma_vals(1) ...
                sigma_vals(2) sigma_vals(2) sigma_vals(2) sigma_vals(2) sigma_vals(2) sigma_vals(2)];
    gamma_row =[gamma_vals(1) gamma_vals(1) gamma_vals(1) gamma_vals(2) gamma_vals(2) gamma_vals(2)...
                gamma_vals(1) gamma_vals(1) gamma_vals(1) gamma_vals(2) gamma_vals(2) gamma_vals(2)];
    alpha_row =[alpha_vals(1) alpha_vals(2) alpha_vals(3) alpha_vals(1) alpha_vals(2) alpha_vals(3) ...
                alpha_vals(1) alpha_vals(2) alpha_vals(3) alpha_vals(1) alpha_vals(2) alpha_vals(3)];

    col_names = ["\sigma", "\gamma", "\alpha", "SP -> SP_star", "SP -> SP_bar", "SP -> LS", "LS_det -> LS"];

    Table = table(sigma_row', gamma_row', alpha_row', sp_sp_star', sp_sp_bar', sp_ls', ls_det_LS', 'VariableNames', col_names);

end  %gen_table