clear; clc; format long;
close all;

ka = 1;          
ln = 1;          
K_A = 0.5;       
K_R = 0.5;       
k_A_on = 0.3;    
k_A_off = 1;     
k_R_off = 0.1;   
k_m_off = 2;     
beta = 10;       
n = 8;           
alpha = 30;      
Fsp0 = 0.2;      
Vpm = 10;        
Vd = 3;          
theta = pi/4;    
k_bar = 10;      

k_s = ka * (k_bar - 1) + 1e-6;  
k_star = ka * k_s / (ka + k_s); 
Fp0 = (1 - Vd/Vpm) * Fsp0;      
t0 = Fp0 / cos(theta);          
C_R0 = 0.5;                     
C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / ...
       (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n)); 
kr_on = (k_R_off + k_m_off * (2 / (1 + exp(beta * t0 / k_star / ln)))) ...
        * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0); 
le0 = ln + t0 / k_star;         
lambda0 = le0 / ln;             
tau1 = ln * cos(theta) / (Vpm - Vd); 
tau2 = Fsp0 / (k_s * Vpm);      


Na = 10;                        
d0 = 1500; % nm                     
[podoconnect, xpod, ypod, Npod] = PodoConnectivity_hexagon(Na, d0);

theta_rot = -30 * pi / 180; 
xpod_rot = xpod * cos(theta_rot) - ypod * sin(theta_rot);
ypod_rot = xpod * sin(theta_rot) + ypod * cos(theta_rot);

center_nodes_idx = find(abs(ypod) < 1e-3);
[spatial_pos_nm, sort_idx] = sort(xpod(center_nodes_idx));
sorted_center_nodes = center_nodes_idx(sort_idx);

mid_idx = ceil(length(sorted_center_nodes) / 2);
kymo_nodes = sorted_center_nodes(mid_idx-2 : mid_idx+2); 
pos_5_um = spatial_pos_nm(mid_idx-2 : mid_idx+2) / 1000; 

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
Dc = d0^2 / 8;
Diff_Op = (Dc / d0^2) * L_mat;

tmax = 60;    
dt = 0.002;   
t_vec = 0:dt:tmax;
Nsteps = length(t_vec);
noise_levels = 0.01;
labels = {'Noise = 0 (Deterministic)', 'Noise = 0.05 (Moderate)', 'Noise = 0.15 (High)'};

time_mask = t_vec >= 59.0;
t_kymo = t_vec(time_mask);

CR_kymo_all = cell(length(noise_levels), 1);
global_min = inf;
global_max = -inf;

disp('Starting 60-second SDE Simulations & Video Generation...');
for i = 1:length(noise_levels)
    amp = noise_levels(i);
    fprintf('\n-> Processing Noise Level: %.2f\n', amp);
    

    [CR_history, ~, ~] = run_sde(Npod, Nsteps, dt, amp, ...
        C_R0, C_A0, lambda0, Diff_Op, kr_on, K_R, K_A, k_R_off, k_m_off, k_A_on, k_A_off, beta, n, alpha, tau1, tau2, k_bar);
    
    skip_frames = round(0.1 / dt); 

    frame_idx = 1 + skip_frames; 

    CR_kymo = CR_history(kymo_nodes, time_mask);
    CR_kymo_all{i} = CR_kymo;
    
    global_min = min(global_min, min(CR_kymo(:)));
    global_max = max(global_max, max(CR_kymo(:)));

    video_filename = sprintf('Podosome_Array_Noise_%.2f.mp4', amp);
    fprintf('   Generating video: %s\n', video_filename);
    
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30; 
    v.Quality = 95;
    open(v);

    fig_movie = figure('Name', 'Video Render', 'Color', 'w', 'Position', [100, 100, 800, 800], 'Visible', 'off');
    axes('Position', [0.05 0.05 0.9 0.9]); 
    hold on;

    h_scatter = scatter(xpod_rot, ypod_rot, 350, CR_history(:, 1), 'filled', 'MarkerEdgeColor', 'none');
    
    colormap(jet); 
    axis equal; axis off;
    
    for k = 1:skip_frames:Nsteps
        h_scatter.CData = CR_history(:, k);
        drawnow;
        writeVideo(v, getframe(fig_movie));
    end
    
    close(v); 
    close(fig_movie);
end
disp('All videos and PDF frames generated successfully.');

function [CR_hist, CA_hist, L_hist] = run_sde(Npod, Nsteps, dt, noise_amp, ...
    C_R0, C_A0, lambda0, Diff_Op, kr_on, K_R, K_A, k_R_off, k_m_off, k_A_on, k_A_off, beta, n, alpha, tau1, tau2, k_bar)
    
    CR = ones(Npod, 1) * C_R0;
    CA = ones(Npod, 1) * C_A0;
    L_var = ones(Npod, 1) * lambda0;
    
    CR_hist = zeros(Npod, Nsteps);
    CA_hist = zeros(Npod, Nsteps);
    L_hist  = zeros(Npod, Nsteps);
    
    CR_hist(:, 1) = CR;
    CA_hist(:, 1) = CA;
    L_hist(:, 1)  = L_var;
    
    sqrt_dt = sqrt(dt);
    

    for k = 1:Nsteps-1

        diff_CR = Diff_Op * CR;
        
        f_lambda = 2 ./ (1 + exp(beta * (L_var - 1)));
        Hill_R = (CR.^n) ./ (K_R^n + CR.^n);
        Hill_A = (CR.^n) ./ (K_A^n + CR.^n);
        
        drift_CR = diff_CR + kr_on * Hill_R .* (1 - CR) - k_R_off * (1 + (k_m_off/k_R_off) * f_lambda) .* CR;
        drift_CA = k_A_on * (1 - CA) .* Hill_A - k_A_off * CA;
        
        num = 1/tau1 + (1/tau2 + alpha * drift_CA) .* (L_var - 1);
        den = 1 + k_bar * exp(-alpha * (CA - C_A0));
        drift_L = 1/tau1 - num ./ den;

        if noise_amp > 0
            CR = CR + drift_CR * dt + noise_amp * CR .* randn(Npod, 1) * sqrt_dt;
        else
            CR = CR + drift_CR * dt;
        end
        CA = CA + drift_CA * dt;
        L_var = L_var + drift_L * dt;
        
        CR = max(CR, 1e-6);
        CA = max(CA, 1e-6);

        CR_hist(:, k+1) = CR;
        CA_hist(:, k+1) = CA;
        L_hist(:, k+1)  = L_var;
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