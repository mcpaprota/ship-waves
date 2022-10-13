function NumModel = numModel(domainLength,domainBreadth,domainDepth,timeIncrement,...
    nEigenvaluesI,nEigenvaluesJ,isBoundaryReflective)
% generates structure array containing numerical model parameters
% Author: Maciej Paprota
NumModel.domainLength = domainLength; % (m)
NumModel.domainBreadth = domainBreadth; % (m)
NumModel.domainDepth = domainDepth; % (m)
NumModel.timeIncrement = timeIncrement; % (m)
NumModel.nEigenvaluesI = nEigenvaluesI;
NumModel.nEigenvaluesJ = nEigenvaluesJ;
NumModel.isBoundaryReflective = isBoundaryReflective; % - 
% 0 - periodic, 1 - reflective boundary conditions
end