function f = gen_plots(C, c, theta, tau_ratio, params, ss, periods)

    % Generate plots replicated from the paper (Figure 1)
    % Input:
    %   C: Consumption with tax adjustment
    %   c: Consumption without tax adjustment
    %   theta: Technology 
    %   tau_ratio: Tax ratio
    %   periods: number of periods to plot.

    % Set up the grid
    t =  0:periods-1;
    % Crop vectors to the correct length
    C = C(1:periods);
    c = c(1:periods);
    theta = theta(1:periods);
    tau_ratio = tau_ratio(1:periods);

    % Unpack steady state values
    C_ss = ss(1);
    tau_ss = ss(2);
    tau_ratio_ss= tau_ss / ( 1 - tau_ss);
    theta_ss = params.theta_bar;

    % Create basic plot
    f = figure;
    hold on
    C_line     = line(t  , 100 * (C - C_ss) / C_ss);
    c_line     = line(t  , 100 * (c - C_ss) / C_ss);
    theta_line = line(t  , 100 * (theta - theta_ss) / theta_ss);
    tau_line   = line(t  , 100 * (tau_ratio - tau_ratio_ss) / tau_ratio_ss);   
    
    % Adjust line properties (functional)
    set(C_line, 'Color', "#D81B60")
    set(c_line, 'Color', "#1E88E5")
    set(theta_line, 'LineStyle', ':', "Color", "#004D40")
    set(tau_line, 'LineStyle', '--', 'Color', '#FFC107')
    % Adjust line properties (aesthetics)
    set(C_line, 'LineWidth', 2)
    set(c_line, 'LineWidth', 2)
    set(theta_line, 'LineWidth', 2)
    set(tau_line, 'LineWidth', 2) 
    % Add labels
    hTitle = title("");
    hXLabel = xlabel("γ = " + params.gamma);
    hYLabel = ylabel('Deviation in percent from steady state');
    % Add legend
    hLegend = legend([C_line, c_line, theta_line, tau_line], ...
    "C with tax adjustment", "C without tax adjustment", ...
    "Productivity", "Tax Ratio", "Location", "NorthEast");

    % Adjust font
    set(gca, 'FontName', 'Helvetica')
    set([hTitle, hXLabel, hYLabel], 'FontName', 'Helvetica')
    set([hLegend, gca], 'FontSize', 10, 'FontName', 'Helvetica')
    set([hXLabel, hYLabel], 'FontSize', 12, 'FontWeight' , 'bold')
    set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    % Adjust axes properties
    set(gca, ...
        'Box', 'off', ...
        'GridColor', [.6     .6 .6], ...
        'TickDir', 'out', ...
        'TickLength', [.02 .02], ...
        'XMinorTick', 'off', ...
        'YMinorTick', 'off', ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'XColor', [.3 .3 .3], ...
        'YColor', [.3 .3 .3], ...
        'YTick', 0:0.2:10, ...
        'XTick', 0:5:periods-1, ...
        'Xlim', [0 periods-1], ...
        'LineWidth', 1)
end  %gen_plots