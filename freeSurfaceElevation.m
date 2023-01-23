function eta = freeSurfaceElevation(eta_hat,kappa_i,kappa_j,x,y,isRefl)
% calculates free-surface elevation based on amplitudes and eigenvalues of cosine expansion
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A Fourier Galerkin method for ship waves. Ocean Engineering (accepted manuscript)
% Input data:
% eta_hat - solution amplitudes (m)
% kappa_i, kappa_j - solution aigenvalues (1/m)
% x, y - spatial coordinates (m)
% Output data:
% eta - free surface elevation (m)
if isRefl
    eta = cos(y'*kappa_j')*eta_hat*cos(kappa_i'*x); % Eq. (45)
else
    eta = real(exp(1i*y'*kappa_j')*eta_hat*exp(1i*kappa_i'*x)); % Eq. (8)
end
