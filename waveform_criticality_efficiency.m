function waveform_criticality_efficiency
    
    tic; clc; clear; format long;
    
    ka = 1; ln = 1;
    K_A = 0.5; K_R = 0.5;
    k_A_on = 0.3; k_A_off = 1;
    k_R_off = 0.1; k_m_off = 2;
    beta = 10; n = 8;
    gamma = k_m_off / k_R_off;
    Fsp0 = 0.2; Vpm = 8; Vd = 3; theta = pi/4;
    alpha = 30; 
    
    tmax = 500; 
    dt = 0.001; 
    
    nx = 100;
    k_bar_scan = logspace(log10(1), log10(800), nx);
    
    Asymmetry = zeros(1, nx);
    Amplitude = zeros(1, nx);
    
    fprintf('Running stiffness scan...\n');
    
    parfor i = 1:nx
        kb = k_bar_scan(i);
        [time, trace, is_osc] = run_single_trace(kb, alpha, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                         k_R_off, k_m_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt);
        if is_osc
            Amplitude(i) = (max(trace) - min(trace)) / 2;
            Asymmetry(i) = calculate_last_cycle_asymmetry(time, trace);
        else
            Amplitude(i) = 0;
            Asymmetry(i) = NaN;
        end
    end
    
    log_k = log10(k_bar_scan);
    dAmp_dk = gradient(Amplitude, log_k);
    Susceptibility = abs(dAmp_dk);
    
    Efficiency = Amplitude .* Susceptibility;
    
    axis_font_size = 14;   
    label_font_size = 16;  
    title_font_size = 18;  
    axis_line_width = 2.0; 

    fig2 = figure('Units', 'pixels', 'Position', [1000, 100, 600, 500], 'Color', 'w');
    
    plot(k_bar_scan, Efficiency, 'g-', 'LineWidth', 2.5, 'Color', [0 0.6 0]); 
    hold on;
    fill([k_bar_scan, fliplr(k_bar_scan)], [Efficiency, zeros(size(Efficiency))], ...
         [0 0.6 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
     
     
    set(gca, 'XScale', 'log');
    ylabel('Efficiency E (A \times \chi)', 'FontWeight', 'bold', 'FontSize', label_font_size);
    xlabel('Stiffness \kappa_{bar}', 'FontWeight', 'bold', 'FontSize', label_font_size);
    title('Communication Efficiency', 'FontSize', title_font_size, 'FontWeight', 'bold');
    
    [max_E, idx_E] = max(Efficiency);
    plot(k_bar_scan(idx_E), max_E, 'gp', 'MarkerSize', 16, 'MarkerFaceColor', 'g');
    text(k_bar_scan(idx_E), max_E*1.08, sprintf('Optimal \\kappa \\approx %.1f', k_bar_scan(idx_E)), ...
         'HorizontalAlignment', 'center', 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', axis_font_size);
    
    grid off; box on;
    set(gca, 'LineWidth', axis_line_width, 'FontSize', axis_font_size, 'TickDir', 'in');
    
    fprintf('Saving Figure 2 (C) as HQ PDF...\n');
    set(fig2, 'PaperPositionMode', 'auto');
    print(fig2, 'Fig_Efficiency_C.pdf', '-dpdf', '-r600', '-painters');
    
    toc;
end

function ratio = calculate_last_cycle_asymmetry(time, trace)
    prominence = (max(trace) - min(trace)) * 0.05; 
    [~, locs_p] = findpeaks(trace, 'MinPeakProminence', prominence);
    [~, locs_v] = findpeaks(-trace, 'MinPeakProminence', prominence);
    
    if length(locs_p) < 2 || isempty(locs_v)
        ratio = NaN; return;
    end
    
    idx_peak_end = locs_p(end);
    valid_valleys = locs_v(locs_v < idx_peak_end);
    if isempty(valid_valleys); ratio = NaN; return; end
    idx_valley = valid_valleys(end);
    
    valid_peaks_start = locs_p(locs_p < idx_valley);
    if isempty(valid_peaks_start); ratio = NaN; return; end
    idx_peak_start = valid_peaks_start(end);
    
    t_p_start = time(idx_peak_start);
    t_valley  = time(idx_valley);
    t_p_end   = time(idx_peak_end);
    
    T_fall = t_valley - t_p_start;
    T_rise = t_p_end - t_valley;
    
    if T_fall <= 0 || T_rise <= 0
        ratio = NaN; 
    else
        ratio = T_rise / T_fall;
    end
end

function [t_out, lambda_trace, is_osc] = run_single_trace(k_bar, alpha, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                    k_R_off, k_m_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt)
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    C_R0 = 0.5;
    Fp0 = (1 - Vd/Vpm) * Fsp0;
    t0 = Fp0 / cos(theta);
    le0 = ln + t0 / k_star;
    lambda0 = le0 / ln;
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
    tau1 = ln * cos(theta) / (Vpm - Vd);
    tau2 = Fsp0 / (k_s * Vpm);
    kr_on_const = (k_R_off + k_m_off * (2 / (1 + exp(beta * (lambda0 - 1))))) ...
                  * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
                  
    N_steps = ceil(tmax / dt);
    rec_start_idx = ceil(100 / dt);
    trace_len = N_steps - rec_start_idx;
    trace = zeros(trace_len, 1);
    
    cr = C_R0 + 0.05; ca = C_A0; lambda = lambda0;
    
    idx = 1;
    for i = 1:N_steps
        f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
        
        drift_cr = kr_on_const * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
        drift_ca = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca; 
        
        num = (1 / tau1 + (1 / tau2 + alpha * drift_ca) * (lambda - 1));
        den = (1 + k_bar * exp(-alpha * (ca - C_A0)));
        drift_lambda = 1 / tau1 - num / den;
        
        cr = cr + drift_cr * dt;
        ca = ca + drift_ca * dt;
        lambda = lambda + drift_lambda * dt;
        
        if i > rec_start_idx
            trace(idx) = lambda;
            idx = idx + 1;
        end
    end
    
    lambda_trace = trace;
    t_out = (rec_start_idx + (1:length(trace)))' * dt;
    
    range_val = max(lambda_trace) - min(lambda_trace);
    
    if range_val < 1e-4
        is_osc = false;
    else
        [pks, ~] = findpeaks(lambda_trace, 'MinPeakProminence', range_val * 0.1);
        [vlys, ~] = findpeaks(-lambda_trace, 'MinPeakProminence', range_val * 0.1);
        if length(pks) >= 2 || length(vlys) >= 2
            is_osc = true;
        else
            is_osc = false;
        end
    end
end