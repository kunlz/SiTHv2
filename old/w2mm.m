
function [mm] = w2mm(LE, Ta)

% latent heat of vaporization [J/kg]
lambda = 2.501e6-2361.*Ta;

% Potential Transpiration and Soil evaporation, mm/day
mm = (LE./lambda).*24.*3600;

end