% ------------------------ %
% Soil moisture Constrains %
% ------------------------ %
function [f_sm, f_sm_s] = swc_stress(wa, soilpar, pET, pftpar)
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

% soil parameters
theta_fc = soilpar(5);
theta_wp = soilpar(7);
% theta_c = soilpar(6);

% Canopy height
CH = pftpar(4);
CH_scalar = CH^0.5;
% scale [1, 25] to [1, 5] 
CH_scalar = 4*((CH_scalar-0.7)/4.3)+1;

a = 0.1;
p = 1./(1+pET) - a.*(1./(1+CH));

theta_wpCH = theta_wp/CH_scalar;

% critical soil mositure for different PFTs
theta_c = (1-p).*(theta_fc-theta_wpCH) + theta_wpCH; 

if theta_c < theta_wpCH
    theta_c = theta_wpCH;
elseif theta_c > theta_fc
    theta_c = theta_fc;
end
% theta_c   = soilpar(6);


if wa <= theta_wpCH
    f_sm = 0;
elseif wa >= theta_c
    f_sm = 1;
else
    f_sm = 1 - ((theta_c-wa)./(theta_c-theta_wpCH)).^CH_scalar; 
end


% constraint for soil evaporation
theta_wp_soil = 0;
if wa <= theta_wp_soil
    f_sm_s = 0;
elseif wa >= theta_fc
    f_sm_s = 1;
else
%     f_sm_s = ((wa - theta_wp)./(theta_fc - theta_wp)).^1;
    f_sm_s = (wa-theta_wp_soil)/(theta_fc - theta_wp_soil); % for soil evaporation only
end

end





% --- old version
% % wc = (theta_c - theta_wp) / (theta_fc - theta_wp);
% if wa <= theta_wp
%     f_sm = 0;
% elseif wa >= theta_c
%     f_sm = 1;
% else
% %     f_sm = 1 - (1 - wa ./ wc).^2;
%     f_sm = 1 - ((theta_c-wa) ./ (theta_c-theta_wp)).^2;
% end
% --- old version

% f_sm_s = wa ./ wc;







