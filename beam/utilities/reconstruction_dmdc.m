function [temporal_recon] = reconstruction_dmdc(dmdc_modes, eigenvals, init_cond, init_amp, t_vec, dt, Abar, Bbar, control,subset)

  % Produce temporal reconstruction at times provided in t_vec
  % Input:
  %  dmdc_modes: DMDc modes to be used in reconstruction
  %  eigenvales: eigenvalues associated with provided DMDc modes
  %  init_cond: Initial condition of system
  %  init_amp: Initial amplitudes
  %  t_vec: Vector of times to evaluate continuous reconstruction at
  %  dt: timestep between snapshots provided to DMDc
  %  Abar: DMDc-learned linear approximation of internal dynamics matrix
  %  Bbar: DMDc-learned linear approximation of exogenous forcing matrix
  %  control: time-history of exogeneous forcing/control
  %  subset: modes to be used in computing the 
  % Output:
  %  temporal_recon: reconstructed solution at times match those in t_vec

  if (nargin == 9)
    s = size(dmdc_modes);
    subset = 1:s(1);
  end

  dmdc_modes = dmdc_modes(:,subset);
  %init_cond = init_cond(subset);

  omega = diag(eigenvals); 
  omega = omega(subset);
  omega = log(omega) ./ dt;
  dmdc_modes_pinv = pinv(dmdc_modes);

  control_matrix = inv(Abar) * Bbar;

  N = length(t_vec);
  for k = 1:N
    disp(strcat(['Time-step ', num2str(k), '/', num2str(N)]));
    t = t_vec(k);
    Nt = floor(t / dt);

    % Construct homogeneous part
    temporal_recon_k = dmdc_modes * diag(exp(omega .* t)) * dmdc_modes_pinv * init_cond;

    % Construct inhomogeneous part
    for j = 1:Nt
      temporal_recon_k = temporal_recon_k + dmdc_modes * diag(exp(omega .* (t-(j-1)*dt))) * dmdc_modes_pinv * control_matrix * control(:,j);
    end

    % Add remaining exogeneous term if time is inbetween snapshots
    if ((t/dt) > Nt)
      temporal_recon_k = temporal_recon_k + ((t - Nt*dt)/dt) .* dmdc_modes * diag(exp(omega .* (t-Nt*dt))) * dmdc_modes_pinv * control_matrix * control(:,Nt);
    end
    temporal_recon(:,k) = temporal_recon_k;
  end

end
