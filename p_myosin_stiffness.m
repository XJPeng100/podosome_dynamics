function p_myosin_stiffness

    tic; clc; clear; format long;
    
    ka = 1; ln = 1; K_A = 0.5; K_R = 0.5;
    k_A_off = 1; k_R_off = 0.1; k_m_off = 2; 
    beta = 10; n = 8; gamma = k_m_off / k_R_off;
    Fsp0 = 0.2; Vpm = 8; Vd = 3; theta = pi/4; alpha = 40;
    
    tmax = 300; dt = 0.01;
    
    nk = 200;
    k_bar_scan = logspace(log10(1), log10(1000), nk);
    
    n_ka = 200;
    ka_on_scan = linspace(0.2, 2, n_ka);
    
    ka_on_specific = [0.3, 0.6, 0.9];
    num_specific = length(ka_on_specific);
    
    fprintf('Running 2D parametric sweep for Heatmap...\n');
    total_iters_2d = n_ka * nk;
    Period_Array_2d = NaN(total_iters_2d, 1);
    
    parfor idx = 1:total_iters_2d
        [i, j] = ind2sub([n_ka, nk], idx);
        current_kA_on = ka_on_scan(i);
        current_kbar = k_bar_scan(j);
        
        [~, trace, is_osc] = run_single_trace(current_kbar, current_kA_on, alpha, ka, ln, K_A, K_R, k_A_off, ...
                                         k_R_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt);
                                         
        if is_osc
            trace_detrend = detrend(trace);
            [~, locs] = findpeaks(trace_detrend, 'MinPeakProminence', 0.05, 'MinPeakDistance', round(0.2/dt));
            if length(locs) >= 3
                periods = diff(locs) * dt;
                Period_Array_2d(idx) = mean(periods); 
            end
        end
    end
    Period_Matrix = reshape(Period_Array_2d, [n_ka, nk]);
    
    fprintf('Running exact 1D sweeps for Line Plots...\n');
    total_iters_1d = num_specific * nk;
    Freq_Array_1d = NaN(total_iters_1d, 1);
    
    parfor idx = 1:total_iters_1d
        [i, j] = ind2sub([num_specific, nk], idx);
        current_kA_on = ka_on_specific(i);
        current_kbar = k_bar_scan(j);
        
        [~, trace, is_osc] = run_single_trace(current_kbar, current_kA_on, alpha, ka, ln, K_A, K_R, k_A_off, ...
                                         k_R_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt);
                                         
        if is_osc
            trace_detrend = detrend(trace);
            [~, locs] = findpeaks(trace_detrend, 'MinPeakProminence', 0.05, 'MinPeakDistance', round(0.2/dt));
            if length(locs) >= 3
                periods = diff(locs) * dt;
                Freq_Array_1d(idx) = 1 / mean(periods); 
            end
        end
    end
    Freq_Specific_Matrix = reshape(Freq_Array_1d, [num_specific, nk]);
    
    figure('Units', 'pixels', 'Position', [100, 100, 900, 450], 'Color', 'w');
    
    subplot(1, 2, 1);
    Freq_Matrix_Heatmap = 1 ./ Period_Matrix;
    Freq_Matrix_Heatmap(isnan(Freq_Matrix_Heatmap)) = 0; 
    
    imagesc(log10(k_bar_scan), ka_on_scan, Freq_Matrix_Heatmap);
    set(gca, 'YDir', 'normal'); 
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Frequency (min^{-1})';
    c.Label.FontWeight = 'bold';
    
    xticks_val = [2, 10, 50, 200, 500];
    xticks(log10(xticks_val));
    xticklabels(string(xticks_val));
    
    xlabel('Substrate Stiffness \kappa_{bar} (log scale)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('p-myosin activity rate (k_A^{on})', 'FontSize', 12, 'FontWeight', 'bold');
    title('(A) Extended Oscillatory Regime', 'FontSize', 13);
    
    subplot(1, 2, 2);
    hold on;
    
    colors = {'b', 'k', 'r'};
    for i = 1:num_specific
        f_line = Freq_Specific_Matrix(i, :);
        valid_idx = ~isnan(f_line); 
        
        plot(k_bar_scan(valid_idx), f_line(valid_idx), '-', 'Color', colors{i}, ...
             'LineWidth', 2, 'DisplayName', sprintf('k_A^{on} = %.1f', ka_on_specific(i)));
    end
    
    set(gca, 'XScale', 'log');
    xlabel('Substrate Stiffness \kappa_{bar}', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Oscillation Frequency (min^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
    title('(B) Effect of p-myosin on Oscillation Frequency', 'FontSize', 13);
    legend('Location', 'northeast');
    grid on; box on;
   
    sgtitle('Supplementary Figure S10: Coupling of Substrate Stiffness and p-myosin', 'FontSize', 14, 'FontWeight', 'bold');
    
    exportgraphics(gcf, 'Fig_S10_p_myosin_coupling.pdf', 'ContentType', 'vector');
    
    toc;
end

function [t_out, lambda_trace, is_osc] = run_single_trace(k_bar, k_A_on, alpha, ka, ln, K_A, K_R, k_A_off, ...
                                    k_R_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt)
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    C_R0 = 0.5; Fp0 = (1 - Vd/Vpm) * Fsp0; t0 = Fp0 / cos(theta);
    le0 = ln + t0 / k_star; lambda0 = le0 / ln;
    
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
    
    tau1 = ln * cos(theta) / (Vpm - Vd); tau2 = Fsp0 / (k_s * Vpm);
    kr_on_const = (k_R_off + k_R_off*gamma * (2 / (1 + exp(beta * (lambda0 - 1))))) * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
    
    N_steps = ceil(tmax / dt);
    cr = C_R0 + 0.05; ca = C_A0; lambda = lambda0;
    rec_start = floor(N_steps * 0.6);
    trace = zeros(N_steps - rec_start, 1);
    idx = 1;
    
    for i = 1:N_steps
        f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
        drift_cr = kr_on_const * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
        drift_ca = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca;
        num = (1 / tau1 + (1 / tau2 + alpha * drift_ca) * (lambda - 1));
        den = (1 + k_bar * exp(-alpha * (ca - C_A0)));
        drift_lambda = 1 / tau1 - num / den;
        
        cr = max(1e-6, cr + drift_cr * dt); 
        ca = max(1e-6, ca + drift_ca * dt); 
        lambda = lambda + drift_lambda * dt;
        
        if i > rec_start; trace(idx) = lambda; idx = idx + 1; end
    end
    lambda_trace = trace; t_out = (0:length(trace)-1) * dt;
    is_osc = range(lambda_trace) > 1e-3;
end