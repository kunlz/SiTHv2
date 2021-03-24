% ------------------------ %
%  Potential ET partition  %
% ------------------------ %
function [pEc,pEs] = potentialET(Rn,G,LAI,Ta,Pa)
% ------ function input -------
% Ta    :  air temperature, C
% Rn    :  average daily net radiation, W/m^2
% Pa    :  atmospheric pressure, kPa
% LAI   :  Leaf area index, 1
% G     :  Soil heat flux, W/m^2
%
% ------ function output ------
% pEc   :  potnetial Transpiration, mm/day
% pEs   :  potnetial Soil evaporation, mm/day
% ------
k = 0.6; % the empirical extinction coefficient set as 0.6
alpha = 1.26; % PT coefficient for water saturated surface
Cp = 1013; % Specific heat (J kg-1 C-1)
eps = 0.622; % e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)

% Radiation located into soil and canopy, seperately
Rns = exp(-k.*LAI).*Rn;
Rnc = Rn-Rns;

% latent heat of vaporization [J/kg]
lambda = 2.501e6-2361.*Ta;

% Saturation vapour pressure at Ta, (kPa)
es = 0.6108.*exp((17.27.*Ta)./(Ta+237.3));

% Slope of saturation vapour pressure curve at Ta (kPa/degC) 
delta = (4098.*es)./((Ta+237.3).^2); 

% psychrometric constant [Pa/C]
gamma = (Cp.*Pa)./(eps*lambda); % Psychrometric constant (kPa/degC)

% Potential Transpiration and Soil evaporation, mm/day
pEc = ((Rnc.*alpha.*delta./(delta+gamma))./lambda).*24.*3600;
pEs = ((Rns-G).*alpha.*delta./(delta+gamma)./lambda).*24.*3600;
end
