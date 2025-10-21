tic
clc
close all
clear
format long;

%% Preprocessing
% Structural parameters
ka = 1;    % stress fiber initial stiffness
ln = 1;    % stress fiber initial length

% Biochemical signal parameters
K_A = 0.5;         % Apparent dissociation constant of actomyosin
K_R = 0.5;         % Apparent dissociation constant of RhoA-GTP
k_A_on = 0.3;      % Self-activation coefficient of actomyosin
k_A_off = 1;       % Self-deactivation coefficient of actomyosin
k_R_off = 0.1;     % Self-deactivation coefficient of RhoA-GTP
k_m_off = 2;       % Mechanical deactivation coefficient of RhoA-GTP
k_m = 10;          % Sensitivity of mechanical feedback
n = 8;             % Hill coefficient
gamma = k_m_off / k_R_off;
beta = 10;
alpha = 80;        % Sensitivity of chemical feedback

Fsp0 = 0.2;        % Initial force, significantly affects phase diagram
Vpm = 8;           % Initial polymerization speed
Vd = 3;            % Depolymerization speed
theta = pi/4;

tau1 = ln * cos(theta) / (Vpm - Vd);

%% Compute Jacobian Matrix
x0 = 1;    % k_bar starts at 1
x1 = 1000; % k_bar ends at 1000

y0 = 1;
y1 = 100;

z0 = 0.2;
z1 = 0.6;

nn = 600;

% Use logspace for k_bar to maintain logarithmic scaling
xx = logspace(log10(x0), log10(x1), nn);

yy = linspace(y0, y1, nn);
zz = linspace(z0, z1, nn);

[x, y, z] = meshgrid(xx, yy, zz);

% Compute k_s from k_bar: k_bar = 1 + k_s / ka => k_s = (k_bar - 1) * ka
k_bar = x;          % x now represents k_bar
k_s = (k_bar - 1) * ka; % Compute k_s from k_bar
k_star = ka .* k_s ./ (ka + k_s);
tau2 = Fsp0 ./ (k_s .* Vpm);

alpha = y;

Fp0 = (1 - Vd / Vpm) * Fsp0;
t0 = Fp0 / cos(theta);
C_R0 = z;

C_A0 = k_A_on .* C_R0.^n ./ (K_A.^n + C_R0.^n) ./ ...
       (k_A_off + k_A_on .* C_R0.^n ./ (K_A.^n + C_R0.^n));

le0 = ln + t0 ./ k_star;
lambda0 = le0 / ln;

f_lambda0 = 2 / (1 + exp(beta * (lambda0 - 1)));

a11 = k_R_off .* (1 + gamma .* f_lambda0) .* (n .* K_R.^n ./ (K_R.^n + C_R0.^n) - 1 ./ (1 - C_R0));
a13 = k_R_off .* (C_R0 .* beta .* gamma .* f_lambda0 .* (1 - 0.5 .* f_lambda0));
a21 = n .* K_A.^n .* k_A_off .* C_A0 ./ (K_A.^n + C_R0.^n) ./ C_R0;
a22 = -k_A_off ./ (1 - C_A0);
a31 = -alpha .* k_bar ./ (1 + k_bar) .* tau2 ./ tau1 .* a21;
a32 = -alpha .* k_bar ./ (1 + k_bar) ./ tau1 .* (1 + tau2 .* a22);
a33 = -1 ./ (1 + k_bar) ./ tau2;

as = -a11;
bs = -a13;
cs = -a21;
ds = -a22;
es = -a31;
fs = -a32;
gs = -a33;

p = (3 * (as .* gs + gs .* ds - bs .* es + as .* ds) - (as + gs + ds).^2) / 3;
q = (27 * (bs .* cs .* fs - bs .* es .* ds + as .* ds .* gs) - 9 * (as + gs + ds) .* (as .* gs + gs .* ds - bs .* es + as .* ds) + 2 * (as + gs + ds).^3) / 27;


Hopf = as .* gs + gs .* ds - bs .* es + as .* ds - (bs .* cs .* fs - bs .* es .* ds + as .* ds .* gs) ./ (as + ds + gs);
aa = Hopf;


Transcritical = bs .* cs .* fs - bs .* es .* ds + as .* ds .* gs;
bb = Transcritical;


Pitchfork = (q / 2).^2 + (p / 3).^3;
cc = Pitchfork;

for i = 1:nn
    for j = 1:nn
        for k = 1:nn
            if (Hopf(i,j,k) > 0) && (Transcritical(i,j,k) > 0)
                cc(i,j,k) = NaN;
            elseif (Hopf(i,j,k) < 0)
                bb(i,j,k) = NaN;
            elseif (Transcritical(i,j,k) < 0)
                aa(i,j,k) = NaN;
            end
        end
    end
end

figure8 = figure('PaperUnits', 'centimeters');
s1 = isosurface(x, y, z, aa, 0);
p1 = patch(s1);
isonormals(x, y, z, aa, p1);
set(p1, 'FaceColor', [53 76 116]/255, 'FaceAlpha', 0.1);
set(p1, 'EdgeColor', 'none');
hold on
s2 = isosurface(x, y, z, bb, 0);
p2 = patch(s2);
isonormals(x, y, z, bb, p2);
set(p2, 'FaceColor', [110 182 211]/255, 'FaceAlpha', 0.1);
set(p2, 'EdgeColor', 'none');
hold on
s3 = isosurface(x, y, z, cc, 0);
p3 = patch(s3);
isonormals(x, y, z, cc, p3);
set(p3, 'FaceColor', [229 144 40]/255, 'FaceAlpha', 0.1);
set(p3, 'EdgeColor', 'none');
view(125, 15)
set(gca, 'XScale', 'log');

box off

% XYZ range
xL = min(xx);
xR = max(xx);
yL = min(yy);
yR = max(yy);
zL = min(zz);
zR = max(zz);

% Vertices for the 3D box
verts = [xL yL zL;
         xR yL zL;
         xR yR zL;
         xL yR zL;
         xL yL zR;
         xR yL zR;
         xR yR zR;
         xL yR zR];

% Faces of the 3D box
faces = [1 2 3 4;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8;
         5 6 7 8];

% Plot the box
patch('Faces', faces, 'Vertices', verts, 'Facecolor', 'b', 'FaceAlpha', 0.02)

xlim([xL, xR]);
ylim([yL, yR]);
zlim([zL, zR]);

grid off

% Export the figure
exportgraphics(figure8, '3D_plot_kbar.pdf', 'ContentType', 'image', 'Resolution', 600)

toc