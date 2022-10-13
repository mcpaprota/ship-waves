function Constants = constants(gravitationalAcceleration,surfaceTensionCoeff,...
    waterDensity)
% generates structure array containing physical constants
% Author: Maciej Paprota
Constants.gravitationalAcceleration = gravitationalAcceleration; % (m/s^2)
Constants.surfaceTensionCoeff = surfaceTensionCoeff; % (N/m^2)
Constants.waterDensity = waterDensity; % (kg/m^3)
end