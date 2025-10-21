function life_cycle
    tic
    clc
    clear;
    close all;
    format long;
    %% Parameters
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

    tmax = 300;     % Total simulation time to reach stable oscillation
    t_char = 3;     % Characteristic time for scaling
    
    % ODE solver options
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'MaxStep', 0.01);
    %% Preprocessing based on parameters
    % Calculate intermediate parameters
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    C_R0 = 0.5;
    Fp0 = (1 - Vd/Vpm) * Fsp0;
    t0 = Fp0 / cos(theta);
    le0 = ln + t0 / k_star;
    lambda0 = le0 / ln;
    
    % Calculate initial C_A0
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / ...
           (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
    
    % Calculate kr_on
    kr_on = (k_R_off + k_m_off * (2 / (1 + exp(beta * (lambda0 - 1))))) ...
            * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
    
    % Set initial conditions for the ODE solver (with a small perturbation)
    initial_conditions = [C_R0 + 0.01; C_A0; lambda0];

    [t, c] = ode15s(@(t, c) dc(t, c, kr_on, C_A0, K_R, K_A, k_A_on, k_A_off, ...
                             k_R_off, gamma, beta, n, alpha, ln, Fsp0, Vpm, Vd, k_star, k_s, theta, k_bar), ...
                    [0, tmax], initial_conditions, opts);
    

    t_scaled = t_char * t;
    

    idx = t_scaled >= 200 & t_scaled <= 212;
    t_filtered = t_scaled(idx) - 200; 
    c_filtered = c(idx, :);

    C_RR = c_filtered(:,1);
    C_AA = c_filtered(:,2);
    lambda = c_filtered(:,3);
    k_t = k_star.*exp(alpha.*(C_AA - C_A0));
    FF = (lambda - 1) .* k_t; % T_v_f
    VV = Vpm./Vd.*(1 - (lambda - 1).*k_t*cos(theta)/Fsp0); % V_p/V_d
    %% Plotting
    % Figure 1: Time evolution curves
    figure1 = figure('Name', 'Time Evolution (k_bar = 3)', 'Position', [100, 100, 600, 800]);
    
    subplot(5,1,1)
    plot(t_filtered, C_RR, '-b', 'linewidth', 2);
    ylabel('C_R');
    title('System Dynamics for k_bar = 3 (Stable Oscillation)');
    xlim([0 12]);
    %grid on;
    subplot(5,1,2)
    plot(t_filtered, C_AA, '-b', 'linewidth', 2);
    ylabel('C_A');
    xlim([0 12]);
    %grid on;
    subplot(5,1,3)
    plot(t_filtered, lambda * 0.15, '-r', 'linewidth', 2);
    ylabel('l_p+l_d');
    xlim([0 12]);
    %grid on;
    subplot(5,1,4)
    plot(t_filtered, FF, '-k', 'linewidth', 2);
    ylabel('T_{vf}');
    xlim([0 12]);
    %grid on;
    subplot(5,1,5)
    plot(t_filtered, VV, '-k', 'linewidth', 2);
    line([0, 12], [Vd/Vd, Vd/Vd], 'LineStyle', '--', 'Color', 'r');
    ylabel('V_p/V_d');
    xlabel('Time');
    xlim([0 12]);
    %grid on;
    
    % Figure 2: Phase portrait (style from model1.m)
    figure2 = figure('Name', 'Phase Portrait (k_bar = 3)', 'PaperUnits', 'centimeters');
    hold on;
    plot(lambda * 0.15, FF, 'b-', 'linewidth', 2);
    
    title('Phase Portrait for k\_bar = 3');
    xlabel('l_p+l_d', 'FontSize', 18);
    ylabel('T_{vf}', 'FontSize', 18);
    set(gca, 'FontSize', 18, 'Box', 'on');
    xlim([0.183 0.197]);
    ylim([0.173 0.181]);
    %grid on;
    hold off;
    
    % --- ADDED CODE FOR PDF EXPORT ---
    disp('Saving figures as high-resolution PDFs...');
    print(figure1, 'Time_Evolution_k_bar_3', '-dpdf', '-vector', '-r600');
    print(figure2, 'Phase_Portrait_k_bar_3', '-dpdf', '-vector', '-r600');
    disp('Finished saving.');
    % ---------------------------------
    
    toc
end
%% Subfunction: ODE dynamics
function dcdt = dc(~, c, kr_on, C_A0, K_R, K_A, k_A_on, k_A_off, ...
                   k_R_off, gamma, beta, n, alpha, ln, Fsp0, Vpm, Vd, k_star, k_s, theta, k_bar)
    cr = c(1); ca = c(2); lambda = c(3);
    
    tau1 = ln * cos(theta) / (Vpm - Vd);
    tau2 = Fsp0 / (k_s * Vpm);
    f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
    
    dcrdt = kr_on * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
    dcadt = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca;
    dlambdadt = 1 / tau1 - (1 / tau1 + (1 / tau2 + alpha * dcadt) * (lambda - 1)) / (1 + k_bar * exp(-alpha * (ca - C_A0)));
    
    dcdt = [dcrdt; dcadt; dlambdadt];
end