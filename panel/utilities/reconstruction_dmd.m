function [temporal_recon] = reconstruction_dmd(dmd_modes, eigenvals, init_cond, init_amp, t_vec, dt, subset)

  % Produce temporal reconstruction at times provided in t_vec
  % Input:
  %  dmd_modes: DMD modes to be used in reconstruction
  %  eigenvales: eigenvalues associated with provided DMD modes
  %  init_cond: Initial condition of system
  %  init_amp: Initial amplitudes
  %  t_vec: Vector of times to evaluate continuous reconstruction at
  %  dt: timestep between snapshots provided to DMD
  %  subset: modes to be used in computing the reconstruction
  % Output:
  %  temporal_recon: reconstructed solution at times match those in t_vec

  if (nargin == 6)
    s = size(dmd_modes);
    %subset = 1:s(1);
    %subset = 1:min(s);
    subset = 1:s(2);
  end

  %dmd_modes = dmd_modes(:,subset);
  %init_cond = init_cond(subset);

  omega = diag(eigenvals); 
  %omega = omega(subset);
  omega = log(omega) ./ dt;
  dmd_modes_pinv = pinv(dmd_modes);
  alphas = dmd_modes_pinv * init_cond;

  N = length(t_vec);
  for k = 1:N
    disp(strcat(['Time-step ', num2str(k), '/', num2str(N)]));
    t = t_vec(k);
    temporal_recon_k = zeros(length(dmd_modes(:,1)), 1);
    Ns = length(subset); Ns
    for ent = 1:Ns
      j = subset(ent);
      temporal_recon_k = temporal_recon_k + dmd_modes(:,j) .* exp(omega(j)*t) .* alphas(j);
    end
    temporal_recon(:,k) = temporal_recon_k;
  end

%  N = length(t_vec);
%  for k = 1:N
%    disp(strcat(['Time-step ', num2str(k), '/', num2str(N)]));
%    t = t_vec(k);
%    temporal_recon_k = dmd_modes * diag(exp(omega .* t)) * dmd_modes_pinv * init_cond;
%    temporal_recon(:,k) = temporal_recon_k;
%  end

end
