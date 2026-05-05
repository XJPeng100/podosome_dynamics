clear; clc; format long;
close all;

ka = 1;          ln = 1;          
K_A = 0.5;       K_R = 0.5;       
k_A_on = 0.3;    k_A_off = 1;     
k_R_off = 0.1;   k_m_off = 2;     
beta = 10;       n = 8;           
alpha = 60;      Fsp0 = 0.2;      
Vpm = 8;         Vd = 3;          
theta = pi/4;    
Dc = 0.5; 
k_bar_scan = [3,   5,   10,  20,  50,  80,  100, 1000]; 
num_scans = length(k_bar_scan);
d0_conditions = [1, 10]; 
mean_freqs_all = zeros(num_scans, 2);
std_freqs_all = zeros(num_scans, 2);
for cond_idx = 1:2
    d0_current = d0_conditions(cond_idx);

    mean_f = zeros(num_scans, 1);
    std_f  = zeros(num_scans, 1);
    
    parfor idx = 1:num_scans
        kb = k_bar_scan(idx);
       
        freqs = simulate_and_get_frequencies(kb, d0_current, Dc, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                           k_R_off, k_m_off, beta, n, alpha, Fsp0, Vpm, Vd, theta);
                                       
        valid_freqs = freqs(~isnan(freqs));
        
        if ~isempty(valid_freqs)
            mean_f(idx) = mean(valid_freqs);
            std_f(idx)  = std(valid_freqs);
        else
            mean_f(idx) = NaN;
            std_f(idx)  = NaN;
        end
    end
    
    mean_freqs_all(:, cond_idx) = mean_f;
    std_freqs_all(:, cond_idx)  = std_f;
end

figure('Name', 'Podosome Array Frequency Statistics', 'Units', 'pixels', 'Position', [100, 150, 1100, 500], 'Color', 'w');
x_pos = 1:num_scans;
max_val = max(max(mean_freqs_all + std_freqs_all));
if isnan(max_val) || max_val <= 0
    max_val = 1; 
end
y_lim_max = max_val * 1.15;

subplot(1, 2, 1);
hold on;
bar(x_pos, mean_freqs_all(:, 1), 0.6, 'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);
errorbar(x_pos, mean_freqs_all(:, 1), std_freqs_all(:, 1), 'k', ...
    'LineStyle', 'none', 'LineWidth', 2.0, 'CapSize', 10);
xticks(x_pos);
xticklabels(string(k_bar_scan));
xlim([0.2, num_scans + 0.8]);
ylim([0, y_lim_max]); 
xlabel('Substrate Stiffness (k_{bar})', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Oscillation Frequency (min^{-1})', 'FontWeight', 'bold', 'FontSize', 14);
title('Small Spacing', 'FontSize', 16);
grid on; box on;
ax1 = gca;
ax1.XGrid = 'off'; ax1.YGrid = 'on';
set(ax1, 'LineWidth', 1.5, 'FontSize', 12);

subplot(1, 2, 2);
hold on;
bar(x_pos, mean_freqs_all(:, 2), 0.6, 'FaceColor', [0.9 0.4 0.3], 'EdgeColor', 'k', 'LineWidth', 1.5);
errorbar(x_pos, mean_freqs_all(:, 2), std_freqs_all(:, 2), 'k', ...
    'LineStyle', 'none', 'LineWidth', 2.0, 'CapSize', 10);
xticks(x_pos);
xticklabels(string(k_bar_scan));
xlim([0.2, num_scans + 0.8]);
ylim([0, y_lim_max]); 
xlabel('Substrate Stiffness (k_{bar})', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Oscillation Frequency (min^{-1})', 'FontWeight', 'bold', 'FontSize', 14);
title('Large Spacing', 'FontSize', 16);
grid on; box on;
ax2 = gca;
ax2.XGrid = 'off'; ax2.YGrid = 'on';
set(ax2, 'LineWidth', 1.5, 'FontSize', 12);

print(gcf, 'Fig_Frequency_Comparison_d0.pdf', '-dpdf', '-r300');

function freqs = simulate_and_get_frequencies(k_bar, d0, Dc, ka, ln, K_A, K_R, k_A_on, k_A_off, ...
                                       k_R_off, k_m_off, beta, n, alpha, Fsp0, Vpm, Vd, theta)
                                   
    k_s = ka * (k_bar - 1) + 1e-6;  
    k_star = ka * k_s / (ka + k_s); 
    Fp0 = (1 - Vd/Vpm) * Fsp0;      
    t0 = Fp0 / cos(theta);          
    C_R0 = 0.5;                     
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n)); 
    kr_on = (k_R_off + k_m_off * (2 / (1 + exp(beta * t0 / k_star / ln)))) * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0); 
    lambda0 = (ln + t0 / k_star) / ln;             
    tau1 = ln * cos(theta) / (Vpm - Vd); 
    tau2 = Fsp0 / (k_s * Vpm);      
    Na = 10;                        
    [podoconnect, ~, ~, Npod] = PodoConnectivity_hexagon(Na, d0);
    BoundPodo = Npod - (Na-1)*6 + 1:Npod; 
    L_mat = sparse(Npod, Npod);
    for i = 1:Npod
        neighbors = podoconnect{i};
        if length(neighbors) == 6
            for j = 1:length(neighbors)
                L_mat(i, neighbors(j)) = 1;
            end
            L_mat(i, i) = -6;
        end
    end
    Diff_Op = (Dc / d0^2) * L_mat; 
    CR = ones(Npod, 1) * C_R0;
    CA = ones(Npod, 1) * C_A0;
    L_var = ones(Npod, 1) * lambda0;
    CR = CR + 0.05 * C_R0 * (rand(Npod, 1) * 2 - 1);
    CR(BoundPodo) = C_R0; 
    tmax = 100; 
    dt = 0.001; 
    Nsteps = ceil(tmax / dt);
    
    save_dt = 0.001;
    save_interval = round(save_dt / dt);
    N_saved = floor(Nsteps / save_interval) + 1;
    
    CR_hist = zeros(Npod, N_saved);
    CR_hist(:, 1) = CR;
    idx_save = 1;
    for k = 1:Nsteps
        f_lambda = 2 ./ (1 + exp(beta * (L_var - 1)));
        Hill_R = (CR.^n) ./ (K_R^n + CR.^n);
        Hill_A = (CR.^n) ./ (K_A^n + CR.^n);
        
        diff_CR = Diff_Op * CR;
        drift_CR = diff_CR + kr_on * Hill_R .* (1 - CR) - k_R_off * (1 + (k_m_off/k_R_off) * f_lambda) .* CR;
        drift_CA = k_A_on * (1 - CA) .* Hill_A - k_A_off * CA;
        
        num = 1/tau1 + (1/tau2 + alpha * drift_CA) .* (L_var - 1);
        den = 1 + k_bar * exp(-alpha * (CA - C_A0));
        drift_L = 1/tau1 - num ./ den;
        
        CR = CR + drift_CR * dt;
        CA = CA + drift_CA * dt;
        L_var = L_var + drift_L * dt;
        
        CR = max(CR, 1e-6);
        CA = max(CA, 1e-6);
        
        CR(BoundPodo) = C_R0;
        CA(BoundPodo) = C_A0;
        L_var(BoundPodo) = lambda0;
        
        if mod(k, save_interval) == 0
            idx_save = idx_save + 1;
            CR_hist(:, idx_save) = CR;
        end
    end
    
    internal_nodes = setdiff(1:Npod, BoundPodo);
    freqs = NaN(length(internal_nodes), 1);
    
    analyze_start_idx = ceil(N_saved / 2);
    
    for i = 1:length(internal_nodes)
        node_id = internal_nodes(i);
        trace = CR_hist(node_id, analyze_start_idx:end);
        
        trace_range = max(trace) - min(trace);
        if trace_range > 1e-3
            prom = trace_range * 0.1; 
            [~, locs] = findpeaks(trace, 'MinPeakProminence', prom);
            
            if length(locs) >= 2
                period = mean(diff(locs)) * save_dt;
                if period > 0
                    freqs(i) = 1 / period; 
                end
            end
        end
    end
end
function [Podoconnect, xpod, ypod, Npod] = PodoConnectivity_hexagon(Na, d0)
    Npod = 3 * Na * (Na - 1) + 1;
    xpod = zeros(Npod, 1);
    ypod = zeros(Npod, 1);
    Podoconnect = cell(Npod, 1);
    Angle_sides = [pi/3, 0, -pi/3, -2*pi/3, pi, 2*pi/3];
    
    for i = 2:Na
        Nfirst = 3 * (i - 1) * (i - 2) + 2;
        xpod(Nfirst) = -d0 * (i - 1);
        ypod(Nfirst) = 0;
        for j = 1:6 * (i - 1) - 1
            Ang = Angle_sides(floor((j - 1) / (i - 1)) + 1);
            xpod(Nfirst + j) = xpod(Nfirst + j - 1) + d0 * cos(Ang);
            ypod(Nfirst + j) = ypod(Nfirst + j - 1) + d0 * sin(Ang);
        end
    end
    
    for i = 1:Npod
        neighbors = [];
        for j = 1:Npod
            if abs((xpod(i) - xpod(j))^2 + (ypod(i) - ypod(j))^2 - d0^2) < 0.1 * d0
                neighbors = [neighbors, j];
            end
        end
        [~, orderY] = sort(ypod(neighbors));
        if (length(neighbors) == 6) && (xpod(neighbors(orderY(1))) - xpod(neighbors(orderY(5))) ~= 0)
            nn = orderY(1);
            orderY(1) = orderY(2);
            orderY(2) = nn;
        end
        Podoconnect{i} = neighbors(orderY);
    end
end