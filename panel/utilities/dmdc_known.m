function [phi, d, x0, alphas, Abar, Bbar, u_tilde, s_tilde, v_tilde] = dmdc_known(data, control, B_known, svd_rank);

  % Apply DMDc with known control map to data-set
  % Input: 
  %  data: full set of data snapshots
  %  control: control snapshots
  %  B_known: the known or well-approximated control matrix
  % Output: 
  %  phi: DMDc modes
  %  d: DMDc eigenvalues 
  %  x0: initial condition
  %  alphas: initial amplitudes
  %  Abar: DMDc-learned linear approxmiation of internal dynamics matrix
  %  Bbar: DMDc-learned linear approximation of exogenous forcing matrix

  Bbar = B_known;

  input = data(:,1:end-1);
  output = data(:,2:end);

  q = size(input); n = q(1); m = q(2); o = min(n,m);
  q = size(control); l = q(1);

  [u_tilde, s_tilde, v_tilde] = svd(input, 'econ');
  if (nargin > 4)   % If svd_rank provided, truncate to low-order model
    u_tilde = u_tilde(:, 1:svd_rank);
    s_tilde = s_tilde(1:svd_rank, 1:svd_rank);
    v_tilde = v_tilde(:, 1:svd_rank);
  end

%   % Check for singular value == 0
%   tol = 1e-10; check = 1; trunc = 1;
%   while ((check > tol) & (trunc < o))
%     disp(trunc); disp(check);
%     if (s_tilde(trunc,trunc) > tol)
%       check = s_tilde(trunc,trunc);
%       trunc = trunc + 1;
%     else
%       check = 1e-15;
%       trunc = trunc - 1;
%     end
%   end
% 
%   % Truncate singular value decomposition to eliminate modes with singular value = 0
%   u_tilde = u_tilde(:,1:trunc); v_tilde = v_tilde(:,1:trunc);
%   s_tilde = s_tilde(1:trunc, 1:trunc);

  % Construct linear approximation of internal dynamics
  Abar = (output - B_known*control)*v_tilde*inv(s_tilde)*u_tilde';
  A_tilde = (u_tilde')*(output - B_known*control)*v_tilde*inv(s_tilde);
  [w,d] = eig(A_tilde);
  d_diag = diag(d);
  for j = 1:length(d_diag)
    if (norm(d_diag(j)) > 1e-12)
      phi_loc = (output - B_known*control)*v_tilde*inv(s_tilde)*w(:,j);
    else
      phi_loc = u_tilde*w(:,j);
    end
    phi(:,j) = phi_loc;
  end

  % Construct initial data and modal amplitudes
  x0 = input(:,1);
  alphas = pinv(phi)*x0;

end

