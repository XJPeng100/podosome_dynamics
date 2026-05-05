function biochemical_noise
    tic
    clc; clear; close all;
    format long;

    %% ================= Task 1: CR Time Series Comparison =================

    noise_levels_curves = [0, 0.02, 0.04]; 
    colors = {'k', 'b', 'r'}; 
    line_styles = {'-', '-', '-'};
    
    fig1 = figure('Name', 'Comparison of CR Dynamics', 'Units', 'inches', 'Position', [1, 1, 8, 5], 'Color', 'w');
    hold on;
    
    fprintf('=== Task 1: Generating Curves ===\n');
    for i = 1:length(noise_levels_curves)
        amp = noise_levels_curves(i);
        fprintf('Simulating for Curve Plot: Noise = %.2f...\n', amp);
        
        [t_plot, cr_plot, ~, ~, ~] = core_simulation(amp, 1); 

        plot(t_plot, cr_plot, 'Color', colors{i}, 'LineStyle', line_styles{i}, ...
             'LineWidth', 2.5, 'DisplayName', ['Noise = ' num2str(amp)]);
    end

    xlabel('Time (scaled)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('C_R Concentration', 'FontSize', 12, 'FontWeight', 'bold');
    title('Effect of Noise on C_R Dynamics', 'FontSize', 14);
    legend('show', 'Location', 'northeast', 'FontSize', 10);
    set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'Box', 'on'); 
    xlim([0, 20]); 
    ylim([0.4, 0.6]);
    hold off;

    fprintf('Saving Figure 1 as PDF...\n');
    set(fig1, 'PaperPositionMode', 'auto'); 
    print(fig1, 'Fig1_TimeEvolution', '-dpdf', '-painters', '-r600');

    %% ================= Task 2: Statistics Bar Charts =================
    noise_levels_stats = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
    
    mean_periods = zeros(size(noise_levels_stats));
    std_periods  = zeros(size(noise_levels_stats));
    cv_values    = zeros(size(noise_levels_stats));
    
    fprintf('\n=== Task 2: Calculating Statistics ===\n');
    n_runs_for_stats = 20; 
    
    for i = 1:length(noise_levels_stats)
        amp = noise_levels_stats(i);
        fprintf('Calculating Stats: Noise = %.2f (Running %d sims)...\n', amp, n_runs_for_stats);

        [~, ~, m_period, s_period, cv_val] = core_simulation(amp, n_runs_for_stats);
        
        mean_periods(i) = m_period;
        std_periods(i)  = s_period;
        cv_values(i)    = cv_val;
    end

    fig2 = figure('Name', 'Statistical Analysis', 'Units', 'inches', 'Position', [1, 1, 7, 8], 'Color', 'w');
    
    subplot(2, 1, 1);
    b1 = bar(noise_levels_stats, mean_periods, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'none');
    hold on;

    errorbar(noise_levels_stats, mean_periods, std_periods, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
    
    ylabel('Mean Period (s)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Period Stability vs. Noise', 'FontSize', 14);
    grid on;
    set(gca, 'LineWidth', 1.2, 'FontSize', 12);
    
    ylim([0, max(mean_periods)*1.2]); 
    
    subplot(2, 1, 2);
    b2 = bar(noise_levels_stats, cv_values, 'FaceColor', [0.8, 0.4, 0.4], 'EdgeColor', 'none');
    
    ylabel('CV of Period', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Noise Amplitude (Multiplicative)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Coefficient of Variation (Period Robustness)', 'FontSize', 14);
    grid on;
    set(gca, 'LineWidth', 1.2, 'FontSize', 12);

    xtips = noise_levels_stats;
    ytips = cv_values;
    labels = string(round(cv_values, 4));
    text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
    
    fprintf('Saving Figure 2 as PDF...\n');
    set(fig2, 'PaperPositionMode', 'auto');
    print(fig2, 'Fig2_Statistics', '-dpdf', '-painters', '-r600');
    
    toc
end

%% ================= Core Simulation Function =================
function [t_out, cr_out, mean_period, std_period, cv_period] = core_simulation(noise_amp, n_sims)
    
    ka = 1; ln = 1;
    K_A = 0.5; K_R = 0.5;
    k_A_on = 0.3; k_A_off = 1;
    k_R_off = 0.1; k_m_off = 2;
    beta = 10; n = 8;
    gamma = k_m_off / k_R_off;
    Fsp0 = 0.2;
    Vpm = 8;
    Vd = 3;
    theta = pi/4;
    alpha = 35; 
    k_bar = 3;  
    
    tmax = 400;     
    t_char = 3;     
    dt = 0.005;     
    N_steps = ceil(tmax / dt);
    
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    C_R0 = 0.5;
    Fp0 = (1 - Vd/Vpm) * Fsp0;
    t0 = Fp0 / cos(theta);
    le0 = ln + t0 / k_star;
    lambda0 = le0 / ln;
    
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / ...
           (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
    
    kr_on_const = (k_R_off + k_m_off * (2 / (1 + exp(beta * (lambda0 - 1))))) ...
            * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
    
    tau1 = ln * cos(theta) / (Vpm - Vd);
    tau2 = Fsp0 / (k_s * Vpm);
    sqrt_dt = sqrt(dt);

    all_periods = [];
    
    for sim = 1:n_sims
        c_series = zeros(N_steps + 1, 3);
        current_c = [C_R0 - 0.05; C_A0; lambda0];
        c_series(1, :) = current_c;
        
        dW = randn(N_steps, 2);
        
        for i = 1:N_steps
            cr = current_c(1);
            ca = current_c(2);
            lambda = current_c(3);
            
            f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
            
            drift_cr = kr_on_const * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
            drift_ca = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca;
            
            num = (1 / tau1 + (1 / tau2 + alpha * drift_ca) * (lambda - 1));
            den = (1 + k_bar * exp(-alpha * (ca - C_A0)));
            drift_lambda = 1 / tau1 - num / den;
            
            cr_new = cr + drift_cr * dt + noise_amp * cr * dW(i, 1) * sqrt_dt;
            ca_new = ca + drift_ca * dt; 
            lambda_new = lambda + drift_lambda * dt;
            
            if cr_new < 1e-6; cr_new = 1e-6; end
            if ca_new < 1e-6; ca_new = 1e-6; end
            
            current_c = [cr_new; ca_new; lambda_new];
            c_series(i+1, :) = current_c;
        end

        t_scaled = t_char * (0:N_steps)' * dt;
        
        cutoff = round(N_steps * 0.5);
        t_seg = t_scaled(cutoff:end);
        cr_seg = c_series(cutoff:end, 1);
        
        estimated_T = 8;
        [~, locs_t] = findpeaks(cr_seg, t_seg, ...
            'MinPeakHeight', max(cr_seg)*0.5, ...
            'MinPeakDistance', estimated_T * 0.6);
        
        periods = diff(locs_t);
        
        if ~isempty(periods)
            median_p = median(periods);
            if std(periods) < 1e-9 
                valid_periods = periods;
            else
                idx_valid = periods > median_p * 0.7 & periods < median_p * 1.3;
                valid_periods = periods(idx_valid);
            end
            all_periods = [all_periods; valid_periods];
        end
        
        if sim == n_sims
            plot_idx = t_scaled >= (t_scaled(end) - 1200);
            t_out = t_scaled(plot_idx);
            t_out = t_out - t_out(1); 
            cr_out = c_series(plot_idx, 1);
        end
    end
    

    if isempty(all_periods)
        mean_period = 0; std_period = 0; cv_period = 0;
    else
        mean_period = mean(all_periods);
        std_period = std(all_periods);
        cv_period = std_period / mean_period;
    end
end