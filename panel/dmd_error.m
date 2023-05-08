
close all;
clear all;
clc;

fs = 30;
ffs = 22;
sz = 36;
neg_epsilon = -1e-10;

% Add utilities path
addpath('./utilities');
addpath('./hashimoto_fem/struct_energies_utils');
addpath('./hashimoto_fem/fluid_function_utils');

% Read in data
data_choice = '13'; M = str2num( strcat([ data_choice(1), '.', data_choice(2:end) ]) );
pressure_fname = strcat(['./data/m', data_choice, '/surface_pressure_history_m', data_choice, '.csv']);
x_fname = strcat(['./data/m', data_choice, '/surface_xloc_history_m', data_choice, '.csv']);
y_fname = strcat(['./data/m', data_choice, '/surface_yloc_history_m', data_choice, '.csv']);
z_fname = strcat(['./data/m', data_choice, '/surface_zloc_history_m', data_choice, '.csv']);
xdot_fname = strcat(['./data/m', data_choice, '/surface_xdot_history_m', data_choice, '.csv']);
ydot_fname = strcat(['./data/m', data_choice, '/surface_ydot_history_m', data_choice, '.csv']);
zdot_fname = strcat(['./data/m', data_choice, '/surface_zdot_history_m', data_choice, '.csv']);

data = dlmread(pressure_fname); data = data';
data_x = dlmread(x_fname); data_x = data_x';
data_y = dlmread(y_fname); data_y = data_y';
data_z = dlmread(z_fname); data_z = data_z';
data_xdot = dlmread(xdot_fname); data_xdot = data_xdot';
data_ydot = dlmread(ydot_fname); data_ydot = data_ydot';
data_zdot = dlmread(zdot_fname); data_zdot = data_zdot';

% Determine how many snapshots to use
snap_level = 1;
snapshot_begin = 21;
snapshot_end = 101;
snapshot_vec = snapshot_begin:snap_level:snapshot_end;
data = data(:, snapshot_vec);

% Define simulation parameters
pc2_dt_nondim = 0.001;
amp = 0.0000508;
omega = 691.15038379;
if (M == 1.1)
  Tref = 242.222;
  pref = 34601.89642944;
  rhoref = 0.49774215;
elseif (M == 1.2)
  Tref = 230.05556;
  pref = 28872.53992512;
  rhoref = 0.43729082;
elseif (M == 1.3)
  Tref = 226.88889;
  pref = 35836.7243159;
  rhoref = 0.55034266;
elseif (M >= 1.4)
  Tref = 223.5;
  pref = 41596.29322233;
  rhoref = 0.64847794;
end
gamma = 1.4;
R = pref / (rhoref * Tref);
aref = sqrt(gamma*R*Tref);
Uref = sqrt(pref/rhoref);
t_scale = 1 / Uref;
L = 1;

pc2_dt = pc2_dt_nondim * t_scale;
N = 24000;
snapshot = 240*snap_level;
dmdc_dt = snapshot*pc2_dt;
max_T = N*pc2_dt;
def_shape_fn = dlmread('./hashimoto_fem/Nx21_Ny41_panel_deform.csv');
data_size = size(data); Np = data_size(1); Nt = data_size(2);
Nx = 21; Ny = 41; % Defined as panel dimensions

disp('temp pause'); pause;

% Import structural quantities
elem_node = importdata('./hashimoto_fem/elements_hires.txt'); elem_node = elem_node.data;
node_elem = dlmread('./hashimoto_fem/node_elem_hires.txt');
nodes = dlmread('./hashimoto_fem/nodes_hires.txt');
structural_node_locs = nodes(:,2:4);
surface_node_data = dlmread('./hashimoto_fem/face_nodes_hires.txt');
struct_faces = dlmread('./hashimoto_fem/surf_faces_trim_hires.txt');

mode1_zdir_deform = dlmread('./hashimoto_fem/mode1_zdir_deform_hires.txt');
def_func = mode1_zdir_deform(:,5);
def_func = def_func ./ max(def_func);

% fix location of structural origin
N_nodes = max(size(structural_node_locs));
shift = [zeros(N_nodes, 2), 0.00101854.*ones(N_nodes,1)];
structural_node_locs = structural_node_locs - shift;
N_face = max(size(surface_node_data));
shift = [zeros(N_face, 2), 0.00101854.*ones(N_face,1)];
surface_node_data(:,2:4) = surface_node_data(:,2:4) - shift;

x_loc = data_x(:,1);
y_loc = data_y(:,1);

x_loc_abs = (x_loc - x_loc(1)); x_loc_abs_surf = reshape(x_loc_abs, Nx, Ny); x_loc_abs_surf = x_loc_abs_surf';
y_loc_abs = (y_loc - y_loc(1)); y_loc_abs_surf = reshape(y_loc_abs, Nx, Ny); y_loc_abs_surf = y_loc_abs_surf';

a = x_loc(Nx) - x_loc(1);
x_plot = (x_loc - x_loc(1)) ./ a; x_plot_surf = reshape(x_plot, Nx, Ny); x_plot_surf = x_plot_surf';
y_plot = (y_loc - y_loc(1)) ./ a; y_plot_surf = reshape(y_plot, Nx, Ny); y_plot_surf = y_plot_surf';

% Compute barycentric coordinates
node_bary = zeros(Ny,Nx,5);
z_loc = data_z(:,1);
for j = 1:Ny
  for i = 1:Nx
    zloc_surf(j,i) = z_loc((j-1)*Nx + i);
  end
end

for j = 1:Ny
  for i = 1:Nx
    disp(strcat(['i = ', num2str(i), ', j = ', num2str(j)]));
    if ((i==1) || (i==Nx))
      node_bary(j,i,:) = zeros(1,5);
    elseif ((j==1) || (j==Ny))
      node_bary(j,i,:) = zeros(1,5);
    else
      loc = [x_loc_abs_surf(j,i), y_loc_abs_surf(j,i), zloc_surf(j,i)];
      [belong_elem, fluid_node_bary] = find_struct_face(loc, surface_node_data, structural_node_locs, node_elem, elem_node, struct_faces, neg_epsilon);
      node_bary(j,i,:) = [fluid_node_bary, belong_elem];
    end
  end
end

% Compute generalized coordinate information
q = (data_z(431,:) .* L) ./ def_shape_fn(21,11);
qdot = (data_zdot(431,:) .* Uref) ./ def_shape_fn(21,11); qdot(1) = amp*omega;
qddot = (-omega^2).*q;

% Construct LPT approximated control matrix
B_lpt = zeros(Np, 2*Np);
for j = 1:Ny
  for i = 1:Nx
    if ((j==1) || (j==Ny) || (i==1) || (i==Nx))
      n = [0 0 1];
      dndq = [0 0 0];
    else
      fluid_node_bary = node_bary(j,i,1:4);
      belong_elem = node_bary(j,i,5);
      loc = [x_loc_abs_surf(j,i), y_loc_abs_surf(j,i), zloc_surf(j,i)];
      belong_elem_nodes = elem_node(belong_elem, :);
      belong_elem_deform = 0*def_func(belong_elem_nodes);
      belong_node_locs = structural_node_locs(belong_elem_nodes, :);
      deformed_nodes = [belong_node_locs(:,1:2), belong_node_locs(:,3) + belong_elem_deform];
      n = determine_normal(deformed_nodes, loc, fluid_node_bary);
      loc_def_func = [zeros(10,2), def_func(belong_elem_nodes)];
      dndq = determine_normal_deriv(deformed_nodes, loc_def_func, loc, fluid_node_bary);
    end

    k = (j-1)*Nx + i;
    B_lpt(k,k) = -rhoref*aref*dot([M*aref,0,0], [dndq(1), dndq(2), dndq(3)]) ...
               + rhoref*aref*amp*omega*dot([0,0,def_shape_fn(j,i)], [dndq(1), dndq(2), dndq(3)]);
    B_lpt(k,k+Np) = rhoref*aref*def_shape_fn(j,i);
  end
end
B_lpt = B_lpt ./ pref;
structure = [ones(Nx*Ny,1)*q(snapshot_vec); ones(Nx*Ny,1)*qdot(snapshot_vec)];

err = data - B_lpt*structure; err = err - 1;

show_data = 0;
if (show_data == 1)
  f_data = figure;
  v_fname = strcat(['./animations/dmd_error_modes/m', num2str(data_choice), '_data.mp4']);
  v_data = VideoWriter(v_fname, 'MPEG-4');
  v_data.FrameRate = 10; v_data.Quality = 12; open(v_data);
  p_max = max(max(data)); p_min = min(min(data));
  e_max = max(max(err)); e_min = min(min(err));
  for j = 1:length(snapshot_vec)
  
    subplot(121);
    data_loc = data(:,j); data_loc = reshape(data_loc, Nx, Ny); data_loc = data_loc';
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), data_loc, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, data(:,j), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; caxis([p_min, p_max]); colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$\vec{p} / \vec{p}_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
    subplot(122);
    err_loc = err(:,j); err_loc = reshape(err_loc, Nx, Ny); err_loc = err_loc';
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), err_loc, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, err(:,j), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; caxis([e_min, e_max]); colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$( \vec{p}  - \vec{p}_{in} - B \cdot \vec{u} ) / \vec{p}_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
    sgtitle(strcat(['Pressure State and Error, $M = ', num2str(M), ', t / \Delta t = ', num2str(j+snapshot_begin-1), ' / ', num2str(snapshot_end), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
    frame = getframe(gcf);
    writeVideo(v_data, frame);
  
    clf;
  
  end
  close(v_data);
end

% temp
show_def = 0;
if (show_def == 1)
  f_def = figure;
  z_max = max(max(data_z)); data_z = data_z ./ z_max
  for j = 1:length(snapshot_vec)
  
    data_z_loc = data_z(:,j); data_z_loc = reshape(data_z_loc, Nx, Ny); data_z_loc = data_z_loc';
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), data_z_loc, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, data_z(:,j), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; caxis([0, 1]); colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('First Natural Mode', 'interpreter', 'latex', 'fontsize', fs);
 
    pause;
    clf;

  end
end



%% Perform DMD on error data
[phi_dmd, d_dmd, x0_dmd, alphas_dmd, Abar_dmd] = dmd(err);

%% Obtain eigenvalues magnitudes
d_dmd_diag = diag(d_dmd);
eig_mag = vecnorm(d_dmd_diag, 2, 2); max_eig_mag = max(eig_mag);
disp(strcat(['Max eigenvalue norm: ', num2str(max_eig_mag)]));

%% Obtain singular values
[u,s,v] = svd(err(:,1:end-1), 'econ');
f_singval = figure;
s_diag = diag(s);
semilogy(1:length(s_diag), s_diag, 'o'); grid on;
xlabel('Diagonal Entry $\#$', 'interpreter', 'latex', 'fontsize', fs);
ylabel('$\log \Sigma_{i,i}$', 'interpreter', 'latex', 'fontsize', fs);
title(strcat(['Singular Values of Error, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
sing_val_fname = strcat(['~/individual_meetings/fellows_08032022/figures/dmd_error_panel/m', data_choice, '_singular_values.png']);
%saveas(f_singval, sing_val_fname);

% Plot eigenvalue sorting
t = size(d_dmd); Ne = min(t);
figure;
subplot(131);
plot(1:Ne, real(diag(d_dmd)), 'ob'); grid on;
ylabel('Real( $\lambda$ )', 'interpreter', 'latex', 'fontsize', fs);
title('Real Part of Eigenvalues', 'interpreter', 'latex', 'fontsize', fs);
subplot(132);
plot(1:Ne, imag(diag(d_dmd)), 'ob'); grid on;
ylabel('Imag( $\lambda$ )', 'interpreter', 'latex', 'fontsize', fs);
title('Imag. Part of Eigenvalues', 'interpreter', 'latex', 'fontsize', fs);
subplot(133);
mag = sqrt( real(diag(d_dmd)).^2 + imag(diag(d_dmd)).^2 );
plot(1:Ne, mag, 'ob'); grid on;
ylabel('Mag', 'interpreter', 'latex', 'fontsize', fs);
title('Mag of Eigenvalues', 'interpreter', 'latex' ,'fontsize', fs);

% Find eigenvalues of Abar
eigA = eig(Abar_dmd);
f_eigs = figure; lh = zeros(1,8); f_eigs.Position = [643 -887 877 674];
lh(1) = plot(eigA, 'o'); grid on; hold on;
sf = sqrt( (tan(omega*pc2_dt*snapshot)^2) / (1 + tan(omega*pc2_dt*snapshot)^2) )
lh(2) = yline(sf, 'r', 'DisplayName', 'First Harmonic'); lh(3) = yline(-1*sf, 'r'); 
sf2 = sqrt( (tan(2*omega*pc2_dt*snapshot)^2) / (1 + tan(2*omega*pc2_dt*snapshot)^2) )
lh(4) = yline(sf2, 'b', 'DisplayName', 'Second Harmonic'); lh(5) = yline(-1.*sf2, 'b');
sf3 = sqrt( (tan(3*omega*pc2_dt*snapshot)^2) / (1 + tan(3*omega*pc2_dt*snapshot)^2) )
lh(6) = yline(sf3, 'g', 'DisplayName', 'Third Harmonic'); lh(7) = yline(-1.*sf3, 'g');
th = 0:0.1:2*pi; xt = cos(th); yt = sin(th);
lh(8) = plot(xt,yt,'--k'); pbaspect([1 1 1]);
ax = gca; ax.FontSize = 16;
xlabel('Re($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
ylabel('Im($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
legend(lh([2,4,6]), 'interpreter', 'latex', 'fontsize', fs, 'location', 'eastoutside');
title_fname = strcat(['Eigenvalues of $\bar{A}$, w/o Accel., $M=', num2str(M), '$, SnapLev = ', num2str(snap_level)]);
title(title_fname, 'interpreter', 'latex', 'fontsize', fs);
eig_fname = strcat(['~/individual_meetings/fellows_08032022/figures/dmd_error_panel/m', data_choice, '_eigs.png']);
%saveas(f_eigs, eig_fname);


%% Plot modes
show_modes = 0;
if (show_modes == 1)
  phi_surf = zeros(Ny, Nx, Ne);
  for j = 1:Ne
    phi_loc = phi_dmd(:,j);
    phi_loc = reshape(phi_loc, Nx, Ny); phi_loc = phi_loc';
    phi_surf(:,:,j) = phi_loc;
  end
  
  f_modes = figure; f_modes.Position = [44 187 1548 655];
  v_fname = strcat(['./animations/dmd_error_modes/m', num2str(data_choice), '.mp4']);
  v_modes = VideoWriter(v_fname, 'MPEG-4');
  v_modes.FrameRate = 1; v_modes.Quality = 12;
  open(v_modes);
  
  for j = 1:Ne
  
    subplot(131);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), real(phi_surf(:,:,j)), 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, real(phi_dmd(:,j)), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
    title(strcat(['Re( $\phi_{', num2str(j), '}$ )']), 'interpreter', 'latex', 'fontsize', fs);
    
    subplot(132);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), imag(phi_surf(:,:,j)), 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, imag(phi_dmd(:,j)), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title(strcat(['Im( $\phi_{', num2str(j), '}$ )']), 'interpreter', 'latex', 'fontsize', fs);
    sg_name = strcat(['DMD-Learned Modes from Error; $\phi_{', num2str(j), '}, M = ', num2str(M), '$']);
    sgtitle(sg_name, 'interpreter', 'latex', 'fontsize', fs);
  
    subplot(133);
    scatter(real(diag(d_dmd)), imag(diag(d_dmd)), 26, 'b'); grid on; hold on;
    scatter(real(d_dmd(j,j)), imag(d_dmd(j,j)), 26, 'red', 'filled');
    th = 0:0.1:2*pi; xt = cos(th); yt = sin(th);
    plot(xt,yt,'--k'); pbaspect([1 1 1]);
    xlabel('Re($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('Im($\lambda$)', 'interpreter', 'latex', 'fontsize', fs);
    eval_tit = strcat(['Re ( $\lambda_{', num2str(j), '}$ ) = ', num2str(real(d_dmd(j,j))), ', Im( $\lambda_{', num2str(j), '}$ ) = ', num2str(imag(d_dmd(j,j)))]);
    title(eval_tit, 'interpreter', 'latex', 'fontsize', fs);

    sgtitle(strcat(['DMD-Learned Modes of Error, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level)]), 'interpreter', 'latex', 'fontsize', fs);
    frame_name = strcat(['./learned_modes/m', data_choice, '/mode', num2str(j), '.png']);
    %saveas(f_modes, frame_name);
  
    frame = getframe(gcf);
    writeVideo(v_modes, frame);

    pause;
    clf;
  
  end
  close(v_modes);
end

% disp('Computing optimal amplitudes...');
% % Temp: spDMD debugging
% subset = 13:14;
% [u,s,v] = svd(err(:,1:end-1), 'econ');
% r = max(size(s));
% F_dmd = (u')*(err(:,2:end))*v*inv(s);
% [y,mu] = eig(F_dmd);
% mu = diag(mu); mu_t = mu(:).^(Nt-2:-1:0);
% v_and = fliplr(mu_t);
% P = ((y')*y).*(conj(v_and*(v_and')));
% q = conj(diag(v_and*v*(s')*y));
% E = ones(r,Np); val = 1;
% for j = 1:r
%   if(sum(j==subset) == 1)
%     E(j,:) = zeros(1,Np);
%   else
%     set = zeros(1,Np);
%     set(val) = 1;
%     E(j,:) = set;
%     val = val + 1;
%   end
% end
% alpha_temp = [eye(r), zeros(r,Np)]*pinv([P, E; E', zeros(Np,Np)])*[q; zeros(Np,1)];
% alphas = alpha_temp(subset);
% pause;

%% Obtain full order continuous reconstruction
t_vec = (1:(snapshot_end-snapshot_begin+1)) .* (snapshot*pc2_dt); t_vec = t_vec - snapshot*pc2_dt;

% Construct fake acceleration information for now
C = zeros(Np, Np);
accel = [ ones(Np,1)*qddot(snapshot_vec) ];

if (M == 1.1)
  continuous_subset = 1:2; % M=1.1
elseif (M == 1.2)
  continuous_subset = 1:2;
  %continuous_subset = 10:11; % M=1.2
elseif (M == 1.3)
  continuous_subset = 1:2;
  %continuous_subset = 11:12; % M=1.3
elseif (M == 1.4)
  continuous_subset = 3:4; % M=1.4
end
[continuous_temporal_recon] = reconstruction_error_model_continuous(phi_dmd, d_dmd, data(:,1), B_lpt, C, snapshot*pc2_dt, err, structure, accel, ones(Np,1), t_vec);
[continuous_temporal_recon_low] = reconstruction_error_model_continuous(phi_dmd, d_dmd, data(:,1), B_lpt, C, snapshot*pc2_dt, err, structure, accel, ones(Np,1), t_vec, continuous_subset);

t_phys = snapshot_vec .* (snapshot*pc2_dt);
t_phys = t_phys .* (M*sqrt(gamma*R*Tref) / a);

show_full_cont_recon = 0;
if (show_full_cont_recon == 1)
  f_full_cont_recon = figure; f_full_cont_recon.Position = [166 216 1346 547];
  v_fname = strcat(['./animations/dmd_recon/m', num2str(data_choice), '_full_cont_recon_spdmd.mp4']);
  %v_fname = strcat(['~/individual_meetings/fellows_09022022/figures/spDMD_panel_full/m', num2str(data_choice), '_full_cont_recon_spDMD.mp4']);
  v_full_cont_recon = VideoWriter(v_fname, 'MPEG-4');
  v_full_cont_recon.FrameRate = 1; v_full_cont_recon.Quality = 12;
  open(v_full_cont_recon);
  
  for j = 1:length(t_phys)
  
    t_loc = t_phys(j);
    loc_cont_recon = continuous_temporal_recon(:,j);
    loc_cont_recon_surf = reshape(loc_cont_recon, Nx, Ny); loc_cont_recon_surf = loc_cont_recon_surf';
    loc_data_surf = reshape(data(:,j), Nx, Ny); loc_data_surf = loc_data_surf';
  
    max_p_recon = max(max(real(loc_cont_recon_surf))); max_p_data = max(max(loc_data_surf)); max_p = max([max_p_recon, max_p_data]);
    min_p_recon = min(min(real(loc_cont_recon_surf))); min_p_data = min(min(loc_data_surf)); min_p = min([min_p_recon, min_p_data]);

    subplot(131);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), real(loc_cont_recon_surf), 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, real(loc_cont_recon), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); caxis([min_p, max_p]); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('Full Recon.', 'interpreter', 'latex', 'fontsize', fs);
  
    subplot(132);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), loc_data_surf, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, data(:,j), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); caxis([min_p, max_p]); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('Original Data', 'interpreter', 'latex', 'fontsize', fs);
  
    diff_surf = real(loc_cont_recon_surf) - loc_data_surf; diff = real(loc_cont_recon) - data(:,j);
    subplot(133);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), diff_surf, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, diff, 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$\%$ Difference', 'interpreter', 'latex', 'fontsize', fs);
  
    sgtitle(strcat(['Full Continuous Reconstruction vs. Original Data, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level), ', $t U_{\infty} / L = ', num2str(t_loc), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
    frame = getframe(gcf);
    writeVideo(v_full_cont_recon, frame);
  
    clf;
  
  end
  close(v_full_cont_recon);
end

show_low_cont_recon = 1;
if (show_low_cont_recon == 1)
  f_low_cont_recon = figure; f_low_cont_recon.Position = [166 216 1346 547];
  v_fname = strcat(['./animations/dmd_recon_low/m', num2str(data_choice), '_low_cont_recon_spdmd.mp4']);
  %v_fname = strcat(['~/individual_meetings/fellows_09022022/figures/spDMD_panel_low/m', num2str(data_choice), '_low_cont_recon_spDMD.mp4']);
  v_low_cont_recon = VideoWriter(v_fname, 'MPEG-4');
  v_low_cont_recon.FrameRate = 1; v_low_cont_recon.Quality = 12;
  open(v_low_cont_recon);
  
  for j = 1:length(t_phys)
  
    t_loc = t_phys(j);
    loc_cont_recon = continuous_temporal_recon_low(:,j);
    loc_cont_recon_surf = reshape(loc_cont_recon, Nx, Ny); loc_cont_recon_surf = loc_cont_recon_surf';
    loc_data_surf = reshape(data(:,j), Nx, Ny); loc_data_surf = loc_data_surf';
  
    subplot(131);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), real(loc_cont_recon_surf), 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, real(loc_cont_recon), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    ylabel('$y / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$p_{R} / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
    subplot(132);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), loc_data_surf, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, data(:,j), 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$p_{CFD} / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
    diff_surf = real(loc_cont_recon_surf) - loc_data_surf; diff = real(loc_cont_recon) - data(:,j);
    subplot(133);
    surf(x_plot_surf, y_plot_surf, ones(Ny,Nx), diff_surf, 'FaceAlpha', 0.33); hold on;
    scatter(x_plot, y_plot, sz, diff, 'filled'); grid on;
    pbaspect([1 2 1]);
    colorbar; colormap('jet'); view(2);
    ax = gca; ax.FontSize = ffs;
    xlabel('$x / a$', 'interpreter', 'latex', 'fontsize', fs);
    title('$\big( p_{R} - p_{CFD} \big) / p_{in}$', 'interpreter', 'latex', 'fontsize', fs);
  
    sgtitle(strcat(['Low-Order Continuous Reconstruction vs. Original Data, $M = ', num2str(M), '$, SnapLev = ', num2str(snap_level), ', $t U_{\infty} / L = ', num2str(t_loc), '$']), 'interpreter', 'latex', 'fontsize', fs);
  
    frame = getframe(gcf);
    writeVideo(v_low_cont_recon, frame);
  
    clf;
  
  end
  close(v_low_cont_recon);
end

%%% Additional functions
function [belong_elem, fluid_node_bary] = find_struct_face(fluid_node_loc, surface_node_locs, structural_node_loc, structural_node_elem, structural_elem_node, structural_faces, neg_epsilon)

  % find clsest structural node
  node_dist = fluid_node_loc' - surface_node_locs(:,2:4)';
  node_dist = vecnorm(node_dist);
  [node_dist, sort_ind] = sort(node_dist);

  done = 0;
  for k = 1:length(node_dist)
    closest_node_ind = surface_node_locs(sort_ind(k), 1);

    % go through faces that share node
    elems_containing_node = structural_node_elem(closest_node_ind, :);
    %disp(elems_containing_node)
    for j = 1:length(elems_containing_node)
      % grab an element that contains closest node
      loc_elem = elems_containing_node(j);
  
      % find surface face that belongs to this element
      face_elem = find(structural_faces(:,1) == loc_elem);
      if isempty(face_elem)
        continue;
      end
  
      % get the nodes that comprise the actual surface face belonging to this element
      face_nodes = structural_faces(face_elem, 2:11);
  
      % Identify which corner nodes are on the surface face
      loc_elem_surf_num = find(face_nodes(1:4) == 0);
      loc_elem_surf = [1 2 3 4]; loc_elem_surf(loc_elem_surf_num) = [];
  
      % Backsolve for area coordinates -- in quadratic space
      coords = bary_coords(fluid_node_loc, loc_elem, structural_elem_node, structural_node_loc, loc_elem_surf_num);
      %disp(coords) 
      if ~((coords(1) < neg_epsilon) || (coords(2) < neg_epsilon) || (coords(3) < neg_epsilon) || (coords(4) < neg_epsilon))
        belong_elem = loc_elem;
        fluid_node_bary = coords;
        done = 1;
        break;
      end
    end

    if (done == 1)
      break;
    end
  end

end

function [coordinates] = bary_coords(fluid_node_loc, loc_struc_elem, elem_node, structural_node_loc, loc_elem_surf_num)

  loc_elem_nodes = elem_node(loc_struc_elem, :);
  loc_elem_nodelocs = structural_node_loc(loc_elem_nodes, :);
  fun = @(L) bary_eqns(L, loc_elem_nodelocs, fluid_node_loc);
  if (loc_elem_surf_num == 1)
    L0 = [0, 1/3, 1/3, 1/3];
  elseif (loc_elem_surf_num == 2)
    L0 = [1/3, 0, 1/3, 1/3];
  elseif (loc_elem_surf_num == 3)
    L0 = [1/3, 1/3, 0, 1/3];
  else
    L0 = [1/3, 1/3, 1/3, 0];
  end
  options = optimoptions('fsolve','Display','off');
  coordinates = fsolve(fun, L0, options);

end

function F = bary_eqns(L, loc_elem_nodelocs, fluid_node_loc)

  loc_elem_x = loc_elem_nodelocs(:,1); loc_elem_y = loc_elem_nodelocs(:,2); loc_elem_z = loc_elem_nodelocs(:,3);

  F(1) = loc_elem_x(1)*(2*L(1)-1)*L(1) + loc_elem_x(2)*(2*L(2)-1)*L(2) + loc_elem_x(3)*(2*L(3)-1)*L(3) + loc_elem_x(4)*(2*L(4)-1)*L(4) + ...
               4*loc_elem_x(5)*L(1)*L(2) + 4*loc_elem_x(6)*L(2)*L(3) + 4*loc_elem_x(7)*L(1)*L(3) + ...
               4*loc_elem_x(8)*L(1)*L(4) + 4*loc_elem_x(9)*L(2)*L(4) + 4*loc_elem_x(10)*L(3)*L(4) - fluid_node_loc(1);

  F(2) = loc_elem_y(1)*(2*L(1)-1)*L(1) + loc_elem_y(2)*(2*L(2)-1)*L(2) + loc_elem_y(3)*(2*L(3)-1)*L(3) + loc_elem_y(4)*(2*L(4)-1)*L(4) + ...
               4*loc_elem_y(5)*L(1)*L(2) + 4*loc_elem_y(6)*L(2)*L(3) + 4*loc_elem_y(7)*L(1)*L(3) + ...
               4*loc_elem_y(8)*L(1)*L(4) + 4*loc_elem_y(9)*L(2)*L(4) + 4*loc_elem_y(10)*L(3)*L(4) - fluid_node_loc(2);

  F(3) = loc_elem_z(1)*(2*L(1)-1)*L(1) + loc_elem_z(2)*(2*L(2)-1)*L(2) + loc_elem_z(3)*(2*L(3)-1)*L(3) + loc_elem_z(4)*(2*L(4)-1)*L(4) + ...
                4*loc_elem_z(5)*L(1)*L(2) + 4*loc_elem_z(6)*L(2)*L(3) + 4*loc_elem_z(7)*L(1)*L(3) + ...
                4*loc_elem_z(8)*L(1)*L(4) + 4*loc_elem_z(9)*L(2)*L(4) + 4*loc_elem_z(10)*L(3)*L(4) - fluid_node_loc(3);

  F(4) = L(1) + L(2) + L(3) + L(4) - 1;


end

function [interp_val] = bary_tet_interp(L, val)

  interp_val = val(1)*(2*L(1)-1)*L(1) + val(2)*(2*L(2)-1)*L(2) + val(3)*(2*L(3)-1)*L(3) + val(4)*(2*L(4)-1)*L(4) + ...
               4*val(5)*L(1)*L(2) + 4*val(6)*L(2)*L(3) + 4*val(7)*L(1)*L(3) + ...
               4*val(8)*L(1)*L(4) + 4*val(9)*L(2)*L(4) + 4*val(10)*L(3)*L(4);

end

