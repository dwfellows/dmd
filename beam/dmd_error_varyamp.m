
close all;
clear all;
clc;

fs = 18;

% Add utilities path
addpath('./utilities');

% Define functions associated with mode shape
c1 = 4.73; L_panel=0.3;
psi = @(x,Lp) -1.*(cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1)-sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp)));
psi_max = max(psi(0:0.00001:L_panel,L_panel));
psi = @(x,Lp) -(1./psi_max).*(cos(c1.*x./Lp) - cosh(c1.*x./Lp) - ((cos(c1) - cosh(c1))/(sin(c1)-sinh(c1))).*(sin(c1.*x./Lp) - sinh(c1.*x./Lp)));

dpsi_dx = @(x,Lp) (-1./psi_max).*(-c1/Lp).*( ((cos(c1)-cosh(c1))/(sin(c1)-sinh(c1))).*(cos(c1.*x./Lp) - cosh(c1.*x./Lp)) + sin(c1.*x./Lp) + sinh(c1.*x./Lp) );
n1 = @(q,x,Lp) -q.*dpsi_dx(x,Lp).*(q.^2 .* (dpsi_dx(x,Lp)).^2 + 1).^(-1/2);
n2 = @(q,x,Lp) ((q.^2) .* (dpsi_dx(x,Lp)).^2 + 1).^(-1/2);
dn1_dq = @(q,x,Lp) -1.*dpsi_dx(x,Lp).*((q.^2).*(dpsi_dx(x,Lp)).^2 + 1).^(-3/2);
dn2_dq = @(q,x,Lp) -q.*((dpsi_dx(x,Lp)).^2).*((q.^2).*((dpsi_dx(x,Lp)).^2) + 1).^(-3/2);

% Read in data
data_choice = 'origamp'; data_choice_m = '3';
M = str2num( strcat([ data_choice_m(1), '.', data_choice_m(2:end) ]) );
pressure_fname = strcat(['./master_data/m', data_choice_m, '_varyamp_hires/', data_choice, '/surface_pressure_history_m', data_choice_m '_', data_choice, '.csv']);
x_fname = strcat(['./master_data/m', data_choice_m, '_varyamp_hires/', data_choice, '/surface_xloc_history_m', data_choice_m, '_', data_choice, '.csv']);
y_fname = strcat(['./master_data/m', data_choice_m, '_varyamp_hires/', data_choice, '/surface_yloc_history_m', data_choice_m, '_', data_choice, '.csv']);
xdot_fname = strcat(['./master_data/m', data_choice_m, '_varyamp_hires/', data_choice, '/surface_xdot_history_m', data_choice_m, '_', data_choice, '.csv']);
ydot_fname = strcat(['./master_data/m', data_choice_m, '_varyamp_hires/', data_choice, '/surface_ydot_history_m', data_choice_m, '_', data_choice, '.csv']);

% Temporary for higher-resolution M=0.85 case
% pressure_fname = './master_data/m085_x241_data/surface_pressure_history_m085_x241.csv';
% x_fname = './master_data/m085_x241_data/surface_xloc_history_m085_x241.csv';
% y_fname = './master_data/m085_x241_data/surface_yloc_history_m085_x241.csv';
% xdot_fname = './master_data/m085_x241_data/surface_xdot_history_m085_x241.csv';
% ydot_fname = './master_data/m085_x241_data/surface_ydot_history_m085_x241.csv';

data = dlmread(pressure_fname); data = data';
data_x = dlmread(x_fname); data_x = data_x';
data_y = dlmread(y_fname); data_y = data_y';
data_xdot = dlmread(xdot_fname); data_xdot = data_xdot';
data_ydot = dlmread(ydot_fname); data_ydot = data_ydot';

nondim_factor = 1;

% Determine how many snapshots to use
snap_level = 1;
snapshot_begin = 1001;
snapshot_end = 5001;
snapshot_vec = snapshot_begin:snap_level:snapshot_end;
data = data(:, snapshot_vec);

% Define simulation parameters
pc2_dt_nondim = 4.5970 * 10^(-5);

if (M< 1)
  amp = 0.000722158;
elseif (M>1)
  amp = 0.00336757;
end

if (strcmp(data_choice,'twiceamp'))
  amp = amp*2;
elseif (strcmp(data_choice,'halfamp'))
  amp = amp / 2;
elseif (strcmp(data_choice,'fifthamp'))
  amp = amp / 5;
elseif (strcmp(data_choice,'tenthamp'))
  amp = amp / 10;
end

omega = 376.99111843;
pref = 97326.6883347;
rhoref = 1.2792;
gamma = 1.44;
Tref = 273;
R = pref / (rhoref * Tref);
aref = sqrt(gamma*R*Tref);
Uref = sqrt(pref/rhoref);
t_scale = 1 / Uref;
L = 1;

pc2_dt = pc2_dt_nondim * t_scale; %pc2_dt = 1.66666e-7;
N = 500000;
snapshot = 100*snap_level;
dmdc_dt = snapshot*pc2_dt;
max_T = N*pc2_dt;
def_shape_fn = dlmread('./old_data/x61_beam_mode1.csv');
%def_shape_fn = dlmread('./master_data/m085_x241_data/x121_beam_mode1.csv'); % Temporary for higher resolution case
data_size = size(data); Nx = data_size(1); Nt = data_size(2);
L_panel = 0.3; x = linspace(0,L_panel, Nx);
q = data_y(round(Nx/2), :) .* L ./ psi(x(round(Nx/2)), L_panel);
qdot = (data_ydot(round(Nx/2), :) .* Uref) ./ psi(x(round(Nx/2)), L_panel);
qddot = -omega^2.*q;

figure;
subplot(311); plot(q); grid on;
subplot(312); plot(qdot); grid on;
subplot(313); plot(qddot); grid on;

% Construct LPT approximated control matrix
B_lpt = zeros(Nx, 2*Nx);
for j = 1:Nx
  B_lpt(j,j) = -rhoref*aref*dot([M*aref,0], [dn1_dq(0,x(j),L_panel), dn2_dq(0,x(j),L_panel)]) ...
%             + rhoref*aref*dot([0,psi(x(j),L_panel)],[dn1_dq(0,x(j),L_panel), dn2_dq(0,x(j),L_panel)]);
              + rhoref*aref*amp*omega*dot([0,psi(x(j),L_panel)],[dn1_dq(0,x(j),L_panel), dn2_dq(0,x(j),L_panel)]);
  B_lpt(j,j+Nx) = rhoref*aref*psi(x(j),L_panel);
end
B_lpt = B_lpt ./ pref;

% Compute LPT predictions
lpt_pred = zeros(Nx,length(snapshot_vec));
for j = 1:length(snapshot_vec)

  loc_snap_num = snapshot_vec(j);
  q_loc = q(loc_snap_num); qdot_loc = qdot(loc_snap_num);
  w_mis = (M*sqrt(gamma*R*Tref)).*-n1(q_loc,x,L_panel);
  w_mot = qdot_loc.*(def_shape_fn).*n2(q_loc,x,L_panel);
  w = w_mis + w_mot;
  pn = pref.*ones(1,Nx) + (rhoref*aref).*w;
  pn = pn ./ pref;
  lpt_pred(:,j) = pn';

end

% Construct acceleration approximation matrix
C = zeros(Nx, Nx);
if (M<1)
  syms eta
  for j = 1:Nx
    A0_loc = (1/pi) * int(log(abs(eta)), -x(j)/L_panel, 1 - x(j)/L_panel);
    C(j,j) = rhoref * A0_loc * psi(x(j),L_panel);
  end
else
  C = zeros(Nx,Nx);
end
C = C ./ pref;
%C = zeros(Nx,Nx); 

structure = [ ones(Nx,1)*q(snapshot_vec); ones(Nx,1)*qdot(snapshot_vec) ];
accel = [ ones(Nx,1)*qddot(snapshot_vec) ];
err = data - B_lpt*structure - C*accel; err = err - 1;
err_noaccel = data - B_lpt*structure; err_noaccel = err_noaccel - 1; % Temp

%% Perform DMD on error data
[phi_dmd, d_dmd, x0_dmd, alphas_dmd, Abar_dmd] = dmd(err);
[phi_dmd_noaccel, d_dmd_noaccel, x0_dmd_noaccel, alphas_dmd_nocall, Abar_dmd_noaccel] = dmd(err_noaccel);

%% Obtain eigenvalue magnitudes
d_dmd_diag = diag(d_dmd);
eig_mag = vecnorm(d_dmd_diag, 2, 2); max_eig_mag = max(eig_mag); format long;
disp(strcat(['Max eigenvalue norm: ', num2str(max_eig_mag)]));

amps = pinv(phi_dmd)*data(:,1);

d_dmd_diag = diag(d_dmd_noaccel);
omega_continuous = log(d_dmd_diag) / dmdc_dt;
[min_val, min_ind] = min(abs(abs(imag(omega_continuous)) - omega));
le_ind = [min_ind, min_ind+1];
continuous_subset = le_ind; discrete_subset = le_ind;
dlmwrite(strcat(['~/research/bem/mode1_modes/m', data_choice_m, '_modes_', data_choice, '.csv']), phi_dmd_noaccel(:,le_ind), 'delimiter', ',', 'precision', 12);

%% Perform discrete reconstruction
t_vec = 1:(snapshot_end - snapshot_begin + 1);
if (strcmp(data_choice,'tenthamp'))
  discrete_subset = 1:2;
elseif (strcmp(data_choice,'fifthamp'))
  discrete_subset = 1:2;
elseif (strcmp(data_choice, 'halfamp'))
  discrete_subset = 5:6;
elseif (strcmp(data_choice, 'origamp'))
  discrete_subset = 3:4;
elseif (strcmp(data_choice, 'twiceamp'))
  discrete_subset = 52:53; % 3001:5001
  %discrete_subset = 48:49; % 4001:5001
end

t_phys = snapshot_vec * (snapshot*pc2_dt);
t_phys = t_phys .* (M*sqrt(gamma*R*Tref) / L_panel);

perform_discrete = 0;
if (perform_discrete)
  [discrete_temporal_recon] = reconstruction_error_model_discrete(phi_dmd, d_dmd, data(:,1), B_lpt, 0.*C, err, structure, accel, ones(Nx,1), t_vec);
  [discrete_temporal_recon_noaccel] = reconstruction_error_model_discrete(phi_dmd_noaccel, d_dmd_noaccel, data(:,1), B_lpt, C, err, structure, accel, ones(Nx,1), t_vec);
  [discrete_temporal_recon_noaccel_low] = reconstruction_error_model_discrete(phi_dmd, d_dmd, data(:,1), B_lpt, 0.*C, err, structure, accel, ones(Nx,1), t_vec, discrete_subset);
  [discrete_temporal_recon_low] = reconstruction_error_model_discrete(phi_dmd, d_dmd, data(:,1), B_lpt, C, err, structure, accel, ones(Nx,1), t_vec, discrete_subset);
  
  t_phys = snapshot_vec * (snapshot*pc2_dt);
  t_phys = t_phys .* (M*sqrt(gamma*R*Tref) / L_panel);
  
  f_recon = figure; f_recon.Position = [34 138 1582 736];
  spatial = linspace(0,1,Nx)';
  spatial = spatial * ones(1,Nt);
  temporal = ones(Nx,1) * t_phys;
  
  recon_sp1 = subplot(141);
  min_p = real(min(min(discrete_temporal_recon))); max_p = real(max(max(discrete_temporal_recon)));
  levels = linspace(min_p, max_p, 100);
  contourf(spatial, temporal, real(discrete_temporal_recon), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Discrete Reconstruction', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_sp2 = subplot(142);
  min_p = min(min(data)); max_p = max(max(data));
  levels = linspace(min_p, max_p, 100);
  contourf(spatial, temporal, data, levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Original Data', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_sp3 = subplot(143);
  diff = real(discrete_temporal_recon - data);
  min_diff = min(min(diff)); max_diff = max(max(diff));
  levels = linspace(min_diff, max_diff, 100);
  contourf(spatial, temporal, diff, levels, 'edgecolor', 'none');
  colorbar; caxis([min_diff, max_diff]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('$\big( \vec{p}_{k}^{DR} - \vec{p}_{k} \big) / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_sp4 = subplot(144);
  min_lpt = min(min(lpt_pred)); max_lpt = max(max(lpt_pred)); levels = linspace(min_lpt, max_lpt, 100);
  contourf(spatial, temporal, lpt_pred, levels, 'edgecolor', 'none');
  colorbar; caxis([min_lpt, max_lpt]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('LPT', 'interpreter', 'latex', 'fontsize', fs);
  
  sgtitle(strcat(['$x-t$ Plots of Reconstruction and Data, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
  
  f_comp = figure; f_comp.Position = [34 138 1582 736];
  diff = real(discrete_temporal_recon - discrete_temporal_recon_noaccel);
  min_p = min(min(diff)); max_p = max(max(diff));
  levels = linspace(min_p, max_p, 100);
  contourf(spatial, temporal, real(diff), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Difference in Reconstruction with Acceleration Term', 'interpreter', 'latex', 'fontsize', fs)
  
  f_recon_low = figure; f_recon_low.Position = [34 138 1582 736];
  
  min_p_full = min(min(real(discrete_temporal_recon))); min_p_low = min(min(real(discrete_temporal_recon_low))); min_p = min([min_p_full, min_p_low]);
  max_p_full = max(max(real(discrete_temporal_recon))); max_p_low = max(max(real(discrete_temporal_recon_low))); max_p = max([max_p_full, max_p_low]);
  levels = linspace(min_p, max_p, 100);
  
  recon_low_sp1 = subplot(141);
  contourf(spatial, temporal, real(discrete_temporal_recon_low), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Low-Order DR', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_sp2 = subplot(142);
  contourf(spatial, temporal, real(discrete_temporal_recon), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Full DR', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_sp3 = subplot(143);
  diff = real(discrete_temporal_recon_low - discrete_temporal_recon);
  min_p = min(min(diff)); max_p = max(max(diff));
  levels = linspace(min_p, max_p, 100);
  contourf(spatial, temporal, diff, levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('$\big( \vec{p}^{DR}_{L} - \vec{p}^{DR}_{F} \big) / \vec{p}_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_sp4 = subplot(144);
  min_lpt = min(min(lpt_pred)); max_lpt = max(max(lpt_pred)); levels = linspace(min_lpt, max_lpt, 100);
  contourf(spatial, temporal, lpt_pred, levels, 'edgecolor', 'none');
  colorbar; caxis([min_lpt, max_lpt]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('LPT', 'interpreter', 'latex', 'fontsize', fs);
  
  
  sgtitle(strcat(['Discrete Reconstructions of Pressure with Acceleration, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
  
  f_recon_low_noaccel = figure; f_recon_low_noaccel.Position = [34 138 1582 736];
  
  min_p_full = min(min(real(discrete_temporal_recon_noaccel))); min_p_low = min(min(real(discrete_temporal_recon_noaccel_low))); min_p = min([min_p_full, min_p_low]);
  max_p_full = max(max(real(discrete_temporal_recon_noaccel))); max_p_low = max(max(real(discrete_temporal_recon_noaccel_low))); max_p = max([max_p_full, max_p_low]);
  levels = linspace(min_p, max_p, 100);
  
  recon_low_noaccel_sp1 = subplot(141);
  contourf(spatial, temporal, real(discrete_temporal_recon_noaccel_low), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Low-Order DR', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_noaccel_sp2 = subplot(142);
  contourf(spatial, temporal, real(discrete_temporal_recon_noaccel), levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('Full DR', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_noaccel_sp3 = subplot(143);
  diff = real(discrete_temporal_recon_noaccel_low - discrete_temporal_recon_noaccel);
  min_p = min(min(diff)); max_p = max(max(diff));
  levels = linspace(min_p, max_p, 100);
  contourf(spatial, temporal, diff, levels, 'edgecolor', 'none');
  colorbar; caxis([min_p, max_p]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('$\big( \vec{p}^{DR}_{L} - \vec{p}^{DR}_{F} \big) / \vec{p}_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
  recon_low_noaccel_sp4 = subplot(144);
  min_lpt = min(min(lpt_pred)); max_lpt = max(max(lpt_pred)); levels = linspace(min_lpt, max_lpt, 100);
  contourf(spatial, temporal, lpt_pred, levels, 'edgecolor', 'none');
  colorbar; caxis([min_lpt, max_lpt]); colormap('jet');
  ax = gca; ax.FontSize = 16;
  xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
  title('LPT', 'interpreter', 'latex', 'fontsize', fs);
  
  
  sgtitle(strcat(['Discrete Reconstructions of Pressure w/o Acceleration, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
end

%% Obtain singular values
[u,s,v] = svd(err(:,1:end-1), 'econ');
f_singval = figure;
s_diag = diag(s);
semilogy(1:length(s_diag), s_diag, 'o'); grid on;
xlabel('Diagonal Entry $\#$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\log \Sigma_{i,i}$', 'interpreter', 'latex', 'fontsize', fs);
title(strcat(['Singular Values of Error, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
sing_val_fname = strcat(['~/individual_meetings/fellows_08012022/figures/corrected_mode1_reconstructions/m', data_choice, '_singular_values.png']);
%saveas(f_singval, sing_val_fname);

% Plot eigenvalue sorting
figure;
subplot(131);
plot(1:Nx, real(diag(d_dmd_noaccel)), 'ob'); grid on;
ylabel('Real( $\lambda$ )', 'interpreter', 'latex', 'fontsize', fs);
title('Real Part of Eigenvalues', 'interpreter', 'latex', 'fontsize', fs);
subplot(132);
plot(1:Nx, imag(diag(d_dmd_noaccel)), 'ob'); grid on;
ylabel('Imag( $\lambda$ )', 'interpreter', 'latex', 'fontsize', fs);
title('Imag. Part of Eigenvalues', 'interpreter', 'latex', 'fontsize', fs);
subplot(133);
mag = sqrt( real(diag(d_dmd_noaccel)).^2 + imag(diag(d_dmd_noaccel)).^2 );
plot(1:Nx, mag, 'ob'); grid on;
ylabel('Mag', 'interpreter', 'latex', 'fontsize', fs);
title('Mag of Eigenvalues', 'interpreter', 'latex' ,'fontsize', fs);

% Find eigenvalues of Abar
eigA = eig(Abar_dmd_noaccel);
f_eigs=figure; f_eigs.Position = [585 -848 761 617];
plot(eigA, 'o'); grid on; hold on;
th = 0:0.1:2*pi; xt = cos(th); yt = sin(th);
plot(xt,yt,'--k'); pbaspect([1 1 1]);
xlabel('Re($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Im($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
title_fname = strcat(['Eigenvalues of $\bar{A}$, w/o Accel., $M=', num2str(M), '$, SnapLev = ', num2str(snap_level)]);
title(title_fname, 'interpreter', 'latex', 'fontsize', fs);
eig_fname = strcat(['~/individual_meetings/fellows_08012022/figures/corrected_mode1_reconstructions/m', data_choice, '_eigs.png']);
%saveas(f_eigs, eig_fname);

f_modes = figure; f_modes.Position = [1 52 1680 895];
for j = 1:Nx
  subplot(12,12,j);
  d_loc = d_dmd_noaccel(j,j); g = real(d_loc); omega = imag(d_loc);
  plot(x ./ L_panel, real(phi_dmd_noaccel(:,j)), 'b'); grid on; hold on;
  plot(x ./ L_panel, imag(phi_dmd_noaccel(:,j)), '--k');
  t_name = strcat(['$Re( \lambda ) = ', num2str(g), ', Im( \lambda ) = ', num2str(omega), '$']);
  title(t_name, 'interpreter', 'latex');
end
title_fname = strcat(['DMD Learned Modes, w/o Accel. $M=', num2str(M), '$, SnapLev = ', num2str(snap_level)]);
sgtitle(title_fname, 'interpreter', 'latex', 'fontsize', fs);
modes_fname = strcat(['~/individual_meetings/fellows_08012022/figures/corrected_mode1_reconstructions/m', data_choice, '_modes.png']);
%saveas(f_modes, modes_fname);

%% Perform continuous reconstruction
t_vec = (1:(snapshot_end - snapshot_begin + 1)) .* (snapshot*pc2_dt); t_vec = t_vec - snapshot*pc2_dt;
[continuous_temporal_recon, alphas_full] = reconstruction_error_model_continuous(phi_dmd_noaccel, d_dmd_noaccel, data(:,1), B_lpt, C, snapshot*pc2_dt, err, structure, accel, ones(Nx,1), t_vec, nondim_factor);
[continuous_temporal_recon_low, alphas] = reconstruction_error_model_continuous(phi_dmd_noaccel, d_dmd_noaccel, data(:,1), B_lpt, C, snapshot*pc2_dt, err, structure, accel, ones(Nx,1), t_vec, nondim_factor, continuous_subset);
dlmwrite(strcat(['~/research/bem/mode1_modes/m', data_choice_m, '_amps_', data_choice, '.csv']), alphas, 'delimiter', ',', 'precision', 12);

amplitude = pinv(phi_dmd_noaccel)*data(:,1);

full_recon = 0;
if (full_recon == 1)
    f_recon_cont = figure; f_recon_cont.Position = [34 138 1582 736];
    spatial = linspace(0,1,Nx)';
    spatial = spatial * ones(1,Nt);
    temporal = ones(Nx,1) * t_phys;
    
    min_p_recon = min(min(real(continuous_temporal_recon))); min_p_orig = min(min(data)); min_p = min([min_p_recon, min_p_orig]);
    max_p_recon = max(max(real(continuous_temporal_recon))); max_p_orig = max(max(data)); max_p = max([max_p_recon, max_p_orig]);
    levels = linspace(min_p, max_p, 100);
    
    recon_sp1 = subplot(141);
    contourf(spatial, temporal, real(continuous_temporal_recon), levels, 'edgecolor', 'none');
    colorbar; caxis([min_p, max_p]); colormap('jet');
    ax = gca; ax.FontSize = 16;
    xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
    title('Continuous Reconstruction', 'interpreter', 'latex', 'fontsize', fs);
    
    recon_sp2 = subplot(142);
    contourf(spatial, temporal, data, levels, 'edgecolor', 'none');
    colorbar; caxis([min_p, max_p]); colormap('jet');
    ax = gca; ax.FontSize = 16;
    xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
    title('Original Data', 'interpreter', 'latex', 'fontsize', fs);
    
    recon_sp3 = subplot(143);
    diff = real(continuous_temporal_recon - data);
    min_diff = min(min(diff)); max_diff = max(max(diff));
    levels = linspace(min_diff, max_diff, 100);
    contourf(spatial, temporal, diff, levels, 'edgecolor', 'none');
    colorbar; caxis([min_diff, max_diff]); colormap('jet');
    ax = gca; ax.FontSize = 16;
    xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
    title('$\big( \vec{p}_{k}^{CR} - \vec{p}_{k} \big) / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
    
    recon_sp4 = subplot(144);
    min_lpt = min(min(lpt_pred)); max_lpt = max(max(lpt_pred)); levels = linspace(min_lpt, max_lpt, 100);
    contourf(spatial, temporal, lpt_pred, levels, 'edgecolor', 'none');
    colorbar; caxis([min_lpt, max_lpt]); colormap('jet');
    ax = gca; ax.FontSize = 16;
    xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
    title('LPT', 'interpreter', 'latex', 'fontsize', fs);
    
    sgtitle(strcat(['$x-t$ Plots of Continuous Reconstruction and Data, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
end

f_recon_cont_low = figure; f_recon_cont_low.Position = [34 138 1582 736];
spatial = linspace(0,1,Nx)';
spatial = spatial * ones(1,Nt);
temporal = ones(Nx,1) * t_phys;

min_p_recon = min(min(real(continuous_temporal_recon_low))); min_p_orig = min(min(data)); min_p = min([min_p_recon, min_p_orig]);
max_p_recon = max(max(real(continuous_temporal_recon_low))); max_p_orig = max(max(data)); max_p = max([max_p_recon, max_p_orig]);
levels = linspace(min_p, max_p, 100);

recon_sp1 = subplot(141);
contourf(spatial, temporal, real(continuous_temporal_recon_low), levels, 'edgecolor', 'none');
colorbar; caxis([min_p, max_p]); colormap('jet');
ax = gca; ax.FontSize = 16;
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$t U_{\infty} / L$', 'interpreter', 'latex', 'fontsize', fs);
title('Low-Order CR', 'interpreter', 'latex', 'fontsize', fs);

recon_sp2 = subplot(142);
contourf(spatial, temporal, data, levels, 'edgecolor', 'none');
colorbar; caxis([min_p, max_p]); colormap('jet');
ax = gca; ax.FontSize = 16;
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
title('Original Data', 'interpreter', 'latex', 'fontsize', fs);

recon_sp3 = subplot(143);
diff = real(continuous_temporal_recon_low) - data;
min_diff = min(min(diff)); max_diff = max(max(diff));
levels = linspace(min_diff, max_diff, 100);
contourf(spatial, temporal, diff, levels, 'edgecolor', 'none');
colorbar; caxis([min_diff, max_diff]); colormap('jet');
ax = gca; ax.FontSize = 16;
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
title('$\big( \vec{p}_{L}^{CR} - \vec{p} \big) / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);

recon_sp4 = subplot(144);
min_lpt = min(min(lpt_pred)); max_lpt = max(max(lpt_pred)); levels = linspace(min_lpt, max_lpt, 100);
contourf(spatial, temporal, lpt_pred, levels, 'edgecolor', 'none');
colorbar; caxis([min_lpt, max_lpt]); colormap('jet');
ax = gca; ax.FontSize = 16;
xlabel('$x / L$', 'interpreter', 'latex', 'fontsize', fs);
title('LPT', 'interpreter', 'latex', 'fontsize', fs);


sgtitle(strcat(['Continuous Reconstructions of Pressure w/o Acceleration, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);

low_cont_fname = strcat(['~/individual_meetings/fellows_08012022/figures/corrected_mode1_reconstructions/m', data_choice, '_continuous_low_recon.png']);
%saveas(f_recon_cont_low, low_cont_fname);


