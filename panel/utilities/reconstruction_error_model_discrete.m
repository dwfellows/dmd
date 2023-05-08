function [discrete_temporal_recon] = reconstruction_error_model_discrete(err_modes, err_eigs, state_init, B, C, err_history, u_history, accel_history, state_mean, t_vec, subset)

  % Produce discrete reconstruction at snapshots provided in t_vec
  % Solves problem: p_{k+1} = A*e_{k} + B*u_{k+1} C*a_{k+1} + p_{mean}
  % Input:
  %   err_modes: Error modes obtained from DMD analysis of error e_{k} = p_{k} - B*u_{k} - C*a_{k} - p_{mean}
  %   err_eigs: Error eigenvalues obtained from DMD analysis of error e_{k} = p_{k} - B*u_{k} - C*a_{k} - p_{mean}
  %   state_init: Initial condition, p_{0}
  %   B: Approximated control matrix associated with external forcing/control u_{k}
  %   e_history: Error history, e_{k}
  %   u_history: Control/forcing history, u_{k}
  %   state_mean: Spatial mean of state
  %   t_vec: Vector of times to evaluate discrete reconstruction at
  %   subset: Modes to be used in computing the reconstruction
  % Output:
  %   discrete_temporal_recon: Reconstructed solution at times matching those in t_vec

  if (nargin == 10)
    s = size(err_modes);
    subset = 1:s(2);
  end

  disp(strcat(['Discrete reconstruction utilizing subset: ', num2str(subset)]));

  max_t = max(t_vec);
  s = size(err_modes); Nx = s(1);
  recon = zeros(Nx, max_t);

  recon(:,1) = state_init;
  A = err_modes(:,subset)*err_eigs(subset,subset)*pinv(err_modes(:,subset));
  disp(state_mean)
  for j = 2:max_t
    recon(:,j) = A*recon(:,j-1) - A*B*u_history(:,j-1) -A*C*accel_history(:,j-1) + (eye(Nx) - A)*state_mean ...
               + B*u_history(:,j) + C*accel_history(:,j);
    %recon(:,j) = A*err_history(:,j-1) + B*u_history(:,j) + state_mean;
  end

  discrete_temporal_recon = recon(:,t_vec);

end
