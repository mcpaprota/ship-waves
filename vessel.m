function Vessel = vessel(hullLength,hullBreadth,hullDraught,...
    hullPositionX,hullPositionY,hullAngleTheta,...
    vesselRoutePositionX,vesselRoutePositionY,vesselRouteAngleTheta)
% generates structure array containing vessel data
% scalars or vectors containing data of M vessel hull components
% centered around coordinate system origin
% Author: Maciej Paprota
Vessel.hullLength = hullLength; % (m)
Vessel.hullBreadth = hullBreadth; % (m)
Vessel.hullDraught = hullDraught; % (m)
Vessel.hullPositionX = hullPositionX; % (m)
Vessel.hullPositionY = hullPositionY; % (m)
Vessel.hullAngleTheta = hullAngleTheta; % (rad)
% vectors containing subsequent data of vessel route
Vessel.routePositionX = vesselRoutePositionX; % (m)
Vessel.routePositionY = vesselRoutePositionY; % (m)
Vessel.routeAngleTheta = vesselRouteAngleTheta; % (rad)
end