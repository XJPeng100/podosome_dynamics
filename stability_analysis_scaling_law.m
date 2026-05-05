function stability_analysis_scaling_law
    tic; clc; clear; format long;
    rng(42);
    %% ---------------- 1. Parameters ----------------
    ka = 1; ln = 1;
    K_A = 0.5; K_R = 0.5;
    k_A_on = 0.3; k_A_off = 1;
    k_R_off = 0.1; k_m_off = 2; 
    beta = 10; n = 8;
    gamma = k_m_off / k_R_off;
    Fsp0 = 0.2;
    Vpm = 8; Vd = 3;
    theta = pi/4;
    
    alpha = 40;   
    tmax = 500; 
    dt = 0.001;
    
    x0 = 3; x1 = 1000; nx = 100;  
    k_bar_scan = logspace(log10(x0), log10(x1), nx);
    
    n_targets = 20; 
    raw_targets = logspace(log10(3), log10(900), n_targets); 
    target_k_bars = unique(round(raw_targets));
    
    n_repeats = 40;         
    noise_amp = 0.018;      

    fprintf('Step 1: Calculating deterministic curve for %d points (Raw Data)...\n', nx);
    
    Freq_det = zeros(1, nx); 
    
    parfor i = 1:nx
        kb = k_bar_scan(i);

        [~, freq] = run_simulation(kb, 0, alpha, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                   k_R_off, k_m_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt);
        Freq_det(i) = freq;
    end

    fprintf('Step 2: Calculating stochastic stats for %d target points...\n', length(target_k_bars));
    
    stats_mean = zeros(1, length(target_k_bars));
    stats_std  = zeros(1, length(target_k_bars));
    
    for j = 1:length(target_k_bars)
        kb_target = target_k_bars(j);
        fprintf('  Processing k_bar = %d ...\n', kb_target);
        
        temp_freqs = zeros(1, n_repeats);
        
        parfor r = 1:n_repeats

            [~, freq] = run_simulation(kb_target, noise_amp, alpha, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                       k_R_off, k_m_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt);
            temp_freqs(r) = freq;
        end
        
        valid_freqs = temp_freqs(temp_freqs > 0);
        
        if isempty(valid_freqs)
            stats_mean(j) = 0;
            stats_std(j) = 0;
        else
            stats_mean(j) = mean(valid_freqs);
            stats_std(j) = std(valid_freqs);
        end
    end

    xxx = 0:0.01:1000;

    theory_mask = xxx >= 1 & xxx <= 10;
    x_theory = xxx(theory_mask);
    y_theory = 0.6 * (ka * (x_theory - 1)).^(-0.3);

    figure('Units', 'pixels', 'Position', [100, 100, 800, 600], 'Color', 'w');
    hold on;
    

    plot(x_theory, y_theory, '--r', 'LineWidth', 2.5, 'DisplayName', 'Theory: 0.6 k_s^{-0.3}');
    
    osc_mask = Freq_det > 0;
    plot(k_bar_scan(osc_mask), Freq_det(osc_mask), '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.3], ...
         'DisplayName', ['Deterministic (\alpha=' num2str(alpha) ')']);
    
    e = errorbar(target_k_bars, stats_mean, stats_std, 'o', ...
             'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
             'LineWidth', 1.5, 'Color', 'b', ...
             'DisplayName', ['Stochastic (N=' num2str(n_repeats) ', \sigma=' num2str(noise_amp) ')']);

    set(gca, 'XScale', 'log');
    xlabel('Linker Stiffness k_{bar}', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Frequency (min^-1)', 'FontSize', 14, 'FontWeight', 'bold');
    title(['Robust Oscillation Frequency Analysis'], 'FontSize', 16);

    %legend('Location', 'northeast', 'FontSize', 12);
    
    xlim([1, 1000]); 
    ylim([0, 0.6]); 
    
    box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 12);
    
    print(gcf, 'Fig3C_Clean.pdf', '-dpdf', '-r600', '-painters');
    toc;
end

function [is_oscillating, frequency] = run_simulation(k_bar, noise_amp, alpha, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                       k_R_off, k_m_off, gamma, beta, n, Fsp0, Vpm, Vd, theta, tmax, dt)
    
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    
    C_R0 = 0.5;
    Fp0 = (1 - Vd/Vpm) * Fsp0;
    t0 = Fp0 / cos(theta);
    le0 = ln + t0 / k_star;
    lambda0 = le0 / ln;
    
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / ...
           (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
           
    tau1 = ln * cos(theta) / (Vpm - Vd);
    tau2 = Fsp0 / (k_s * Vpm);
    
    kr_on_const = (k_R_off + k_m_off * (2 / (1 + exp(beta * (lambda0 - 1))))) ...
                  * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
    N_steps = ceil(tmax / dt);
    
    cr = C_R0 + 0.01; 
    ca = C_A0;
    lambda = lambda0;
    
    record_start_idx = floor(N_steps * 0.5); 
    signal_length = N_steps - record_start_idx;
    lambda_signal = zeros(signal_length, 1);
    
    sqrt_dt = sqrt(dt);
    if noise_amp > 0
        dW = randn(N_steps, 2); 
    end
    
    idx_sig = 1;
    for i = 1:N_steps
        f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
        
        drift_cr = kr_on_const * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
        drift_ca = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca;
        
        numerator = (1 / tau1 + (1 / tau2 + alpha * drift_ca) * (lambda - 1));
        denominator = (1 + k_bar * exp(-alpha * (ca - C_A0)));
        drift_lambda = 1 / tau1 - numerator / denominator;
        
        if noise_amp > 0

            cr = cr + drift_cr * dt + noise_amp * cr * dW(i, 1) * sqrt_dt;
            ca = ca + drift_ca * dt;
        else
            cr = cr + drift_cr * dt;
            ca = ca + drift_ca * dt;
        end
        lambda = lambda + drift_lambda * dt;
        
        if cr < 1e-6; cr = 1e-6; end
        if ca < 1e-6; ca = 1e-6; end
        
        if i > record_start_idx
            lambda_signal(idx_sig) = lambda;
            idx_sig = idx_sig + 1;
        end
    end
    
    sig_detrend = detrend(lambda_signal);
    
    range_val = max(sig_detrend) - min(sig_detrend);
    if range_val < 1e-5
        frequency = 0; is_oscillating = false; return;
    end

    smooth_window = round(0.5 / dt); 
    sig_smooth = smoothdata(sig_detrend, 'gaussian', smooth_window);
    
    prominence_thresh = 0.2 * range_val; 
    
    [~, locs] = findpeaks(sig_smooth, ...
                          'MinPeakProminence', prominence_thresh, ...
                          'MinPeakDistance', round(1.0 / dt));
    
    if length(locs) >= 3
        periods = diff(locs) * dt;
        med_T = median(periods);
        valid_periods = periods(periods > 0.6 * med_T & periods < 1.4 * med_T);
        
        if ~isempty(valid_periods)
            avg_period = mean(valid_periods);
            frequency = 1 / avg_period;
            is_oscillating = true;
        else
            frequency = 1 / med_T; 
            is_oscillating = true;
        end
    else

        [~, locs_retry] = findpeaks(sig_smooth, ...
                          'MinPeakProminence', 0.1 * range_val, ... 
                          'MinPeakDistance', round(1.0 / dt));
        
        if length(locs_retry) >= 3
             periods = diff(locs_retry) * dt;
             frequency = 1 / median(periods);
             is_oscillating = true;
        else
            frequency = 0;
            is_oscillating = false;
        end
    end
end