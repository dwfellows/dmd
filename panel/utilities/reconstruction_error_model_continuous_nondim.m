function [temporal_recon] = reconstruction_error_model_continuous_nondim(err_modes, err_eigs, state_init, B, dt, err_history, u_history, state_mean, t_vec, gamma, theta_max, M, subset)

  % Produce continuous reconstruction at snapshots provided in t_vec
  % Solves problem: p_{k+1} = A*e_{k} + B*u_{k+1} + C*a_{k+1} + p_{mean}

  nondim_qty = gamma*theta_max*(M^2 / sqrt(abs(M^2 - 1)));

  if (nargin == 12)
    s = size(err_modes);
    subset = 1:s(2);
  end

  max_t = size(t_vec); max_t = max_t(2);
  s = size(err_modes); Nx = s(1);
  recon = zeros(Nx, max_t);

  % Compute optimal amplitudes utilizing spDMD constraints
  [u,s,v] = svd(err_history(:,1:end-1), 'econ');
  r = max(size(s));
  F_dmd = (u')*(err_history(:,2:end))*v*inv(s);
  [y,mu] = eig(F_dmd);
  mu = diag(mu); mu_t = mu(:).^(max_t-2:-1:0);
  v_and = fliplr(mu_t);
  P = ((y')*y).*(conj(v_and*(v_and')));
  q = conj(diag(v_and*v*(s')*y));
  E = ones(r,Nx); val = 1;
  for j = 1:r
    if(sum(j==subset) == 1)
      E(j,:) = zeros(1,Nx);
    else
      set = zeros(1,Nx);
      set(val) = 1;
      E(j,:) = set;
      val = val + 1;
    end
  end
  alpha_temp = [eye(r), zeros(r,Nx)]*pinv([P, E; E', zeros(Nx,Nx)])*[q; zeros(Nx,1)];
  alphas = alpha_temp(subset);

  disp(alphas);
  %dlmwrite('alphas_m11.csv', alphas, 'delimiter', ',', 'precision', 12);

  A = err_modes*err_eigs*pinv(err_modes);
  cm = B;

  err_modes = err_modes(:,subset);
  omega = diag(err_eigs); omega = omega(subset);
  omega = log(omega) ./ dt;
  err_modes_pinv = pinv(err_modes);

  for k = 1:max_t
    t = t_vec(k);
    Nt = floor(t / dt);

    % Homogeneous part
    %temporal_recon_k = err_modes * diag(exp(omega.*t)) * err_modes_pinv * (state_init - B*u_history(:,1) - state_mean);  % Utilize initial snapshot for amplitudes
    temporal_recon_k = err_modes * diag(exp(omega.*t)) * alphas; % Utilize spDMD for amplitudes


    % Add on remaining parts
    temporal_recon_k = temporal_recon_k .* nondim_qty;
    temporal_recon_k = temporal_recon_k + B*u_history(:,Nt+1);
    %temporal_recon_k = temporal_recon_k + C*accel_history(:,Nt+1);
    temporal_recon_k = temporal_recon_k + state_mean;
    temporal_recon(:,k) = temporal_recon_k;

    %pause;

  end


end






