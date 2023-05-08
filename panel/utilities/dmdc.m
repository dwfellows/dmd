function [phi, d, x0, alphas, Abar, Bbar, u_tilde, s_tilde, v_tilde] = dmdc(data, control);

  % Apply DMDc to data-set
  % Input: 
  %  data: full set of data snapshots
  %  control: control snapshots
  % Output: 
  %  phi: DMDc modes
  %  d: DMDc eigenvalues 
  %  x0: initial condition
  %  alphas: initial amplitudes
  %  Abar: DMDc-learned linear approxmiation of internal dynamics matrix
  %  Bbar: DMDc-learned linear approximation of exogenous forcing matrix

  input = data(:,1:end-1);
  output = data(:,2:end);

  q = size(input); n = q(1);
  q = size(control); l = q(1);

  x0 = input(:,1);
  Omega_mat = [input; control];
  [u_tilde, s_tilde, v_tilde] = svd(Omega_mat, 'econ');
  Gbar = output * v_tilde * inv(s_tilde) * u_tilde';
  u_tilde1 = u_tilde(1:n,:); u_tilde2 = u_tilde(n+1:end,:);
  Abar = output * v_tilde * inv(s_tilde) * u_tilde1';
  Bbar = output * v_tilde * inv(s_tilde) * u_tilde2';
  [u_hat, s_hat, v_hat] = svd(output,'econ');
  A_tilde = u_hat' * Abar * u_hat;
  B_tilde = u_hat' * Bbar;
  [w,d] = eig(A_tilde);
  phi = output * v_tilde * inv(s_tilde) * u_tilde1' * u_hat * w;
  alphas = pinv(phi)*x0;

end

