function coherence_resonance
    tic; clc; clear; close all;
    
    rng(42, 'twister'); 

    noise_level = 0.02;     
    tmax_plot = 10;         

    [t_det, cr_det] = core_simulation(0, tmax_plot);
    
    [t_stoch, cr_stoch] = core_simulation(noise_level, tmax_plot);
    
    fig = figure('Name', 'Comparison_Dynamics', 'Units', 'inches', 'Position', [2, 2, 8, 5], 'Color', 'w');
    hold on; 
    
    plot(t_stoch, cr_stoch, 'Color', [0.7, 0.7, 1.0], 'LineWidth', 2.0, 'DisplayName', 'Stochastic Path');
    plot(t_det, cr_det, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Deterministic Mean');
    
    xlabel('Time (min)', 'FontSize', 12);
    ylabel('C_R Concentration', 'FontSize', 12);
    title('Reproducible Stochastic Dynamics', 'FontSize', 14);
    
    box on;        
    grid on;                
    set(gca, 'FontSize', 13, 'LineWidth', 1.2); 
    
    xlim([0, max(t_stoch)]);
    ylim([0.46, 0.58]);
    
    legend('Location', 'northeast');
    hold off;
    print(fig, 'Comparison', '-dpdf', '-painters', '-r600');

end


function [t_out, cr_out] = core_simulation(noise_amp, tmax_val)
 
    ka=1; ln=1; K_A=0.5; K_R=0.5; k_A_on=0.3; k_A_off=1; k_R_off=0.1; k_m_off=2;
    beta=10; n=8; gamma=k_m_off/k_R_off; Fsp0=0.2; Vpm=8; Vd=3; theta=pi/4;
    alpha=35; k_bar=3; t_char=3; dt=0.005; 
    
    N_steps = ceil(tmax_val / dt);
    k_s = ka * (k_bar - 1) + 1e-6;
    k_star = ka * k_s / (ka + k_s);
    C_R0 = 0.52; 
    Fp0 = (1 - Vd/Vpm) * Fsp0;
    lambda0 = (ln + (Fp0 / cos(theta)) / k_star) / ln;
    C_A0 = k_A_on * C_R0^n / (K_A^n + C_R0^n) / (k_A_off + k_A_on * C_R0^n / (K_A^n + C_R0^n));
    kr_on_const = (k_R_off + k_m_off * (2 / (1 + exp(beta * (lambda0 - 1))))) * (K_R^n + C_R0^n) / C_R0^(n-1) / (1 - C_R0);
    
    tau1 = ln * cos(theta) / (Vpm - Vd);
    tau2 = Fsp0 / (k_s * Vpm);
    

    cr = C_R0 + 0.00; ca = C_A0; lambda = lambda0;
    cr_out = zeros(N_steps + 1, 1);
    cr_out(1) = cr;

    for i = 1:N_steps
        f_lambda = 2 / (1 + exp(beta * (lambda - 1)));
        d_cr = kr_on_const * cr^n / (K_R^n + cr^n) * (1 - cr) - k_R_off * (1 + gamma * f_lambda) * cr;
        d_ca = k_A_on * (1 - ca) * cr^n / (K_A^n + cr^n) - k_A_off * ca;
        d_lambda = 1/tau1 - (1/tau1 + (1/tau2 + alpha * d_ca) * (lambda - 1)) / (1 + k_bar * exp(-alpha * (ca - C_A0)));

        cr = cr + d_cr * dt + noise_amp * cr * randn * sqrt(dt);
        ca = ca + d_ca * dt;
        lambda = lambda + d_lambda * dt;

        if cr < 1e-6; cr = 1e-6; end
        cr_out(i+1) = cr;
    end
    t_out = t_char * (0:N_steps)' * dt;
end