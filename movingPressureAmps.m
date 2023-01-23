function p_hat = movingPressureAmps(kappa_i,kappa_j,l,b,L,B,D,theta,X,Y,X_m0,Y_m0,theta_m0,isRefl)
% moving pressure function spectral amplitudes
% Author: Maciej Paprota
% Reference: M. Paprota. 2023. A Fourier Galerkin method for ship waves. Ocean Engineering (accepted manuscript)
% Input data:
% L - ship hull m-th component length (m)
% B - ship hull m-th component breadth (m)
% D - ship hull m-th component draught (m)
% theta - history of ship rotation angle with respect to x-axis (rad)
% X_m0, Y_m0 - ship hull m-th component center coordinates
% theta_m0 - ship hull m-th component rotation
% when centered around origin of the coordinate system (m)
% X, Y - ship hull system center coordinates (m)
% isRefl - lateral boundary indicator (0 - periodic, 1 - reflective)
p_hat = kappa_i*0;
for m=1:length(X_m0)
    alpha = D(m); 
    sigma_x = L(m)/sqrt(8*log(10)); % Eq. (23)
    sigma_y = B(m)/sqrt(8*log(10)); % Eq. (24)
    beta_x = cos(theta+theta_m0(m)).^2/2/sigma_x^2+sin(theta+theta_m0(m)).^2/2/sigma_y^2; % Eq. (28)
    beta_xy = sin(2*(theta+theta_m0(m)))/4*(1/sigma_x^2-1/sigma_y^2); % Eq. (29)
    beta_y = sin(theta+theta_m0(m)).^2/2/sigma_x^2+cos(theta+theta_m0(m)).^2/2/sigma_y^2; % Eq. (30)
    X_m = X_m0(m)*cos(theta)-Y_m0(m)*sin(theta)+X;
    Y_m = X_m0(m)*sin(theta)+Y_m0(m)*cos(theta)+Y;
    if isRefl
    p_hat = p_hat + 4*pi*alpha/l/b/sqrt(beta_x*beta_y-beta_xy^2)*...
        exp((beta_y*kappa_i.^2+beta_x*kappa_j.^2)/4/(beta_xy^2-beta_x*beta_y)).*...
        (cosh(beta_xy*kappa_i.*kappa_j/2/(beta_xy^2-beta_x*beta_y)).*...
        cos(kappa_i*X_m).*cos(kappa_j*Y_m)+ ...
        sinh(beta_xy*kappa_i.*kappa_j/2/(beta_xy^2-beta_x*beta_y)).*...
        sin(kappa_i*X_m).*sin(kappa_j*Y_m)); % Eqs. (51)
    else
    p_hat = p_hat + pi*alpha/l/b/sqrt(beta_x*beta_y-beta_xy^2)*...
        exp((beta_y*kappa_i.^2-2*beta_xy*kappa_i.*kappa_j+...
        beta_x*kappa_j.^2)/4/(beta_xy^2-beta_x*beta_y)).*...
        (exp(-1i*kappa_i*X_m-1i*kappa_j*Y_m)); % Eqs. (32)
    end
end
if isRefl % Eq. (52)
    p_hat(1,:) = p_hat(1,:)/2;
    p_hat(:,1) = p_hat(:,1)/2;
end
end
