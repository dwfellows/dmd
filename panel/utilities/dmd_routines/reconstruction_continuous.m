function [temporal_recon, alphas] = reconstruction_continuous(modes, eigs, state_init, dt, history, state_mean, t_vec, subset)

  % Inputs:
    % modes: DMD modes
    % eigs: DMD eigenvalues
    % state_init: initial snapshot
    % dt: time between DMD snapshots
    % history: snapshot history
    % state mean: steady state / average state; enter zeros if mean not subtracted off
    % t_vec: vector of times to compute reconstruction
    % subset: vector of modes to use from subset
  % Outputs
    % temporal_recon: reconstruction of data using DMD modes at provided t_vec
    % alphas: DMD modal amplitudes used in reconstruction
    
  % Produce continuous reconstruction at snapshots provided in t_vec
  % Solves problem: p_{k+1} = A*e_{k} + p_{mean}

  if (nargin == 7)
    s = size(modes);
    subset = 1:s(2);
  end

  max_t = size(t_vec); max_t = max_t(2);
  s = size(modes); Nx = s(1);
  recon = zeros(Nx, max_t);

  % % Compute optimal amplitudes utilizing spDMD constraints
  % [u,s,v] = svd(err_history(:,1:end-1), 'econ');
  % r = max(size(s));
  % F_dmd = (u')*(err_history(:,2:end))*v*inv(s);
  % [y,mu] = eig(F_dmd);
  % mu = diag(mu); mu_t = mu(:).^(max_t-2:-1:0);
  % v_and = fliplr(mu_t);
  % P = ((y')*y).*(conj(v_and*(v_and')));
  % q = conj(diag(v_and*v*(s')*y));
  % E = ones(r,Nx); val = 1;
  % for j = 1:r
  %   if(sum(j==subset) == 1)
  %     E(j,:) = zeros(1,Nx);
  %   else
  %     set = zeros(1,Nx);
  %     set(val) = 1;
  %     E(j,:) = set;
  %     val = val + 1;
  %   end
  % end
  % alpha_temp = [eye(r), zeros(r,Nx)]*pinv([P, E; E', zeros(Nx,Nx)])*[q; zeros(Nx,1)];
  % alphas = alpha_temp(subset);

  alphas = pinv(err_modes(:,subset))*err_history(:,1);

  A = modes*eigs*pinv(modes);

  modes = modes(:,subset);
  omega = diag(eigs); omega = omega(subset);
  omega = log(omega) ./ dt;
  modes_pinv = pinv(modes);

  for k = 1:max_t
    t = t_vec(k);
    Nt = floor(t / dt);

    % Homogeneous part
    temporal_recon_k = modes * diag(exp(omega.*t)) * alphas; % Utilize projection amplitudes

    % Add on remaining parts
    temporal_recon_k = temporal_recon_k + state_mean;
    temporal_recon(:,k) = temporal_recon_k;

  end


end






