function p = movingPressure(x,y,Vessel,n)
% moving pressure function
% Author: Maciej Paprota
% Input data:
L = Vessel.hullLength;
B = Vessel.hullBreadth;
D = Vessel.hullDraught;
X_m0 = Vessel.hullPositionX;
Y_m0 = Vessel.hullPositionY;
theta_m0 = Vessel.hullAngleTheta;
X = Vessel.routePositionX(n);
Y = Vessel.routePositionY(n);
theta = Vessel.routeAngleTheta(n);
x = ones(length(y),1)*x;
y = y'*ones(1,length(x(1,:)));
p = x*0;
for m=1:length(X_m0)
    alpha = D(m); 
    sigma_x = L(m)/sqrt(8*log(10)); % Eq. (xx)
    sigma_y = B(m)/sqrt(8*log(10)); % Eq. (xx)
    beta_x = cos(theta+theta_m0(m)).^2/2/sigma_x^2+sin(theta+theta_m0(m)).^2/2/sigma_y^2; % Eq. (xx)
    beta_xy = sin(2*(theta+theta_m0(m)))/4*(1/sigma_x^2-1/sigma_y^2); % Eq. (xx)
    beta_y = sin(theta+theta_m0(m)).^2/2/sigma_x^2+cos(theta+theta_m0(m)).^2/2/sigma_y^2; % Eq. (xx)
    X_m = X_m0(m)*cos(theta)-Y_m0(m)*sin(theta)+X;
    Y_m = X_m0(m)*sin(theta)+Y_m0(m)*cos(theta)+Y;
    p = p + alpha*exp(-(beta_x*(x-X_m).^2+2*beta_xy*(x-X_m).*(y-Y_m)+...
        beta_y*(y-Y_m).^2)); % Eqs. (xx) and (xx)
end
end