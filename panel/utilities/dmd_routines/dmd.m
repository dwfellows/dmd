function [phi, d, x0, alphas, Abar] = dmd(data);

  % Apply DMD to data-set
  % Input: 
  %  data: data set snapshots
  % Output: 
  %  phi: DMD modes
  %  d: DMD eigenvalues 
  %  x0: initial condition
  %  alphas: initial amplitudes
  %  Abar: DMD-learned linear approxmiation of internal dynamics matrix

  input = data(:,1:end-1);
  output = data(:,2:end);

  q = size(input); n = q(1);

  x0 = input(:,1);
  [u,s,v] = svd(input, 'econ');
  
  % Check for singular value == 0
%  tol = 1e-10; check = 1; trunc = 1;
%  while (check > tol)
%    disp(trunc); disp(check);
%    if (s(trunc,trunc) > tol)
%      check = s(trunc,trunc);
%      trunc = trunc + 1;
%    else
%      check = 1e-15;
%      trunc = trunc - 1;
%    end
%  end
%
%  u = u(:,1:trunc); v = v(:,1:trunc);
%  s = s(1:trunc, 1:trunc);

  sinv = inv(s);
  Abar = output * v * sinv * u';
  Atilde = (u') * Abar * u;
  [w,d] = eig(Atilde);
  phi = output * v * sinv * w;
  alphas = pinv(phi) * x0;
  % phi_hat = u*w;  % original DMD modes -- see Proctor DMDc paper

end

