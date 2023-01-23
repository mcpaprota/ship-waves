function [kappa, kappa_i, kappa_j, eta_hat, phi_hat, p_hat] = ...
    shipWaves_lm(NumModel,Vessel,Constants)
% Fourier Galerkin solution to gravity-capillary linear ship waves
% low memory version of shipWaves - only a final result is presented
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A Fourier Galerkin method for ship waves
% Input data:
l = NumModel.domainLength;
b = NumModel.domainBreadth;
d = NumModel.domainDepth;
dt = NumModel.timeIncrement;
I = NumModel.nEigenvaluesI;
J = NumModel.nEigenvaluesJ;
isRefl = NumModel.isBoundaryReflective;
L = Vessel.hullLength;
B = Vessel.hullBreadth;
D = Vessel.hullDraught;
X_m0 = Vessel.hullPositionX;
Y_m0 = Vessel.hullPositionY;
theta_m0 = Vessel.hullAngleTheta;
X = Vessel.routePositionX;
Y = Vessel.routePositionY;
theta = Vessel.routeAngleTheta;
g = Constants.gravitationalAcceleration;
gamma = Constants.surfaceTensionCoeff;
rho = Constants.waterDensity;

% Output data:
% kappa_i, kappa_j, kappa - solution eigenvalues (1/m)
% eta_hat - free-surface amplitudes (m)
% phi_hat - velocity potential amplitudes (m^2/s)

% solution eigenvalues
if isRefl
    kappa_i = ones(J,1)*(0:I-1)*pi/l; % Eq. (47)
    kappa_j = (0:J-1)'*ones(1,I)*pi/b; % Eq. (48)
else
    kappa_i = ones(J,1)*(-I/2:I/2-1)*2*pi/l; % Eq. (9)
    kappa_j = (-J/2:J/2-1)'*ones(1,I)*2*pi/b; % Eq. (10)
end
kappa = sqrt(kappa_i.^2+kappa_j.^2); % Eq. (11)
kappatanhkappad = kappa.*tanh(kappa*d);
omega = sqrt((g+gamma*kappa.*kappa/rho)...
    .*kappatanhkappad); % Eq. (38)
% moving pressure term
N = length(X);
p_hat_p = movingPressureAmps(kappa_i,kappa_j,l,b,L,B,D,...
    theta(1),X(1),Y(1),X_m0,Y_m0,theta_m0,isRefl)*g./...
    (g+gamma*kappa.^2/rho);
phi_hat_p = zeros(J,I); % initial velocity potential coefficients - Eq. (17)
eta_hat_p = -p_hat_p; % initial free-surface elevation coefficients - Eq. (18)
% Modified Euler implicit scheme loop:
for n=1:N-1
    p_hat = movingPressureAmps(kappa_i,kappa_j,l,b,L,B,D,...
        theta(n+1),X(n+1),Y(n+1),X_m0,Y_m0,theta_m0,isRefl)*g./...
        (g+gamma*kappa.^2/rho);
    eta_hat = (eta_hat_p-(dt*omega/2).^2.*(eta_hat_p+...
        p_hat+p_hat_p)+dt*phi_hat_p.*kappatanhkappad)./...
        (1+(dt*omega/2).^2); % Eq. (39)
    phi_hat = ((1-(dt*omega/2).^2).*phi_hat_p-...
        dt/2.*(g+gamma*kappa.*kappa/rho).*...
        (2*eta_hat_p+p_hat+p_hat_p))./(1+(dt*omega/2).^2); % Eq. (40);
    eta_hat_p = eta_hat; phi_hat_p = phi_hat; p_hat_p = p_hat;
    disp(['Calculation progress: ' num2str((n+1)/N*100,'%6.2f') '%'])
end
kappa_i = kappa_i(1,:);
kappa_j = kappa_j(:,1);
end
