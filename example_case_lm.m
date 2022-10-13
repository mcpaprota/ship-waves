% Numerical example - low memory version
% circular motion of a catamaran ship
% gravity finite-water waves
% Author: Maciej Paprota
% Reference: M. Paprota. 2022. A Fourier Galerkin method for ship waves
clear, clc, set(0,'defaulttextinterpreter','latex')
% Initialization:
g = 9.8145; % gravitational acceleration (m^2/s)
gamma = 0.074; % surface tension coefficient (N/m^2)
rho = 1000; % water density (kg/m^3)
d = 10; % fluid domain depth (m)
colormap(parula)
Fr_L = 1; % hull Froude number
L = [1 1]; B = [0.3 0.3]; D = [0.5 0.5]; % hull size
X_m0 = [0 0]; Y_m0 = [-1 1]; theta_m0 = [0 0]; % initial hull configuration
isRefl = 1; % define boundaries periodic/reflective 0/1
l = 20; % fluid domain length (m)
b = 20; % fluid domain breadth (m)
I = 100; % number of eigenvalues
J = 100; % number of eigenvalues
% circular motion along circumference
U = sqrt(g*max(L))*Fr_L; % ship target velocity
T = 2*pi/sqrt(2*pi/max(L)*g*tanh(2*pi/max(L)*d)); % ship wave period
dt = T/16; % time increment (s)
NumModel = numModel(l,b,d,dt,I,J,isRefl); % initialize model parameters
r = 5; % circular path radius (m)
simTime = 7*pi/8*r/U; % nearly full circle simulation time
accTime = 0.1*simTime; % (s) acceleration time to target velocity
t = (0:dt:simTime); % time vector (s)
a = U/accTime; % acceleration to target velocity (m^2/s)
sigma = a*t/r; % angular speed in accelerated motion
sigma(t>accTime) = U/r+0*t(t>accTime); % angular speed at target velocity
X = r*sin(sigma.*t)+l/2; % ship route x-coordinate
Y = r*cos(sigma.*t)+b/2; % ship route y-coordinate
theta = -sigma.*t; % rotation of a ship hull according to velocity
Constants = constants(g,gamma,rho); % initialize constants
Vessel = vessel(L,B,D,X_m0,Y_m0,theta_m0,X,Y,theta); % initialize vessel
[kappa, kappa_i, kappa_j, eta_hat, phi_hat, p_hat] = ...
    shipWaves_lm(NumModel,Vessel,Constants); % calculating solution coeffs
% fluid domain
dx = l/1000; x = (0:dx:l);
dy = b/1000; y = (0:dy:b);
% Plotting output:
eta = freeSurfaceElevation(eta_hat,kappa_i,kappa_j,x,y,isRefl);
p = movingPressure(x,y,Vessel,length(X));
p(p<0.01*max(max(p))) = NaN;
graph1 = surf(x,y,eta);
graph1.EdgeColor = 'none';
graph1.FaceAlpha = 0.6;
hold on
graph2 = surf(x,y,-p);
graph2.EdgeColor = 'none';
graph2.FaceColor = 'red';
clim([-1 1]), axis([0 l 0 b -d/2 d/2])
set(gca,'TickLabelInterpreter','latex')
xlabel('$x$'), ylabel('$y$'), zlabel('$\eta$')
set(gca, 'Layer','top')
hold off