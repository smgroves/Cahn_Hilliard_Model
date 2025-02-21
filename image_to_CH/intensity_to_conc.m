indir= cd;

init_file = sprintf("%s/input/MCF10A_T6I8_chromosome_phi_IC_256.csv",indir);
% init_file = sprintf("%s/input/figure_5_chromosome_phi_IC.csv",indir);
phi0 = readmatrix(init_file);
psi0 = (phi0 + 1)/2;
%psi0 is 256X256 matrix with 0.5 for values with 64<x<192  and 0 otherwise
% psi0 = zeros(256,256);
% psi0(64:191,:) = 0.3;
domain = 6.6;
% C_avg = 3.968;
C_avg = 2.757;
% C_avg = 1;
% domain = 3.2;
V_expected = 3.2 * 1.6;  % Original expected volume (for Cavg) in µm^3
V_image = domain*domain/2;  % Actual image volume in µm^3
% Step 1: Compute the corrected average concentration
C_new = C_avg * (V_expected / V_image);

A_image = 3.2 * 1.6;  % Physical area in µm^2
thickness = 1;  % Assumed thickness in µm

% Step 1: Compute volume of image and volume per pixel
V_image = A_image * thickness;  % Total volume in µm^3
num_pixels = numel(psi0);
V_pixel = 2*V_image / num_pixels;  % µm^3 per pixel

% Step 2: Compute sum of intensities
S = sum(psi0(:));

% Step 3: Compute concentration in µmol per pixel
% note V_image and V_pixel are in µm^3; converting to L by 1e-15 gets canceled out when calculating C_uM.
C_pixel = psi0 * (C_new * V_image / S);

% Step 4: Convert to µM (µmol/L)
C_uM = C_pixel / V_pixel;

disp(['Total moles (µmol): ', num2str(sum(C_pixel(:)))]);
disp(['Mean concentration (µM): ', num2str(mean(C_uM(:)))]);

h =heatmap(C_uM);
h.GridVisible = 'off';

writematrix(C_uM,sprintf("%s_C_um_Borealin.csv",init_file));

