% ------------------------ %
% Soil moisture Constrains %
% ------------------------ %
function [S_plant,S_soil] = swc_stress(wa,soilpar)
% -------- function input -------
% wa      : The antecedent soil water content expressed
%           as a function of the WHC in that layer
% soilpar : Soil parameters according to Soil type
% -------- function output -------
% S_plant : Soil moisture stress to plant transpiration
% S_soil  : Soil moisture stress to soil evaporation
% --------
% Stress function for plant transpiration and soil evaporation:
%                  (theta_c-theta_wp)
% wc =  ------------------------------    (about 0.4; Choudhury & Digirolamo, 1998)
%                  (theta_fc-theta_wp)
% where theta_c : the critical soil water content at which plant stress start
% Stree Function (Martens et al., 2017)
%                      theta_c-theta
% S_plant   =   1- (-------------------)^2   =   1-(1-w/wc)^2
%                     theta_c-theta_wp
%
%                    theta_c-theta
% S_soil    =   1- -----------------    =   1-(1-w/wc)=w/wc
%                    theta_c-theta_r
% -------------------------------------------------------------------------

wc = soilpar(8)./soilpar(3);

S_plant_x = 1-(1-wa./wc).^2;
S_soil_x = wa./wc;

S_plant = ones(size(wa));
S_soil = ones(size(wa));

S_plant(wa<wc) = S_plant_x(wa<wc);
S_soil(wa<wc) = S_soil_x(wa<wc);

end