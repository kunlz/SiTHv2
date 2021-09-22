% ----------------------------------------------------------------------- %
%                                                                         %
%                            ==============                               %
%                              SiTHM model                                %
%                            ==============                               %
%                                                                         %
% Written by : Kun Zhang & Gaofeng Zhu;  Lanzhou University, China        %
% Original version : 5/9/2018                                             %
% Last modified : 16/4/2021                                               %
% Questions to zhangkun322@outlook.com                                    %
% ----------------------------------------------------------------------- %


function [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, Pnet, IWS, Vmax] = SiTHM(Rn, Ta, Topt,...
    Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp, optpara) 
% Main function
% -------------------------------------------------------------------------
% Model inout  ::  1  - Rn      -- Net Radiation, W/m-2
%              ::  2  - Ta      -- Near surface air temperature, C
%              ::  2  - Topt    -- Optimal growth temperature for plant, C
%              ::  3  - Pe      -- Precipetition, mm day-1
%              ::  4  - Pa      -- Near surface air pressure, kPa
%              ::  5  - s_VOD   -- Constrains from VOD,[0,1]
%              ::  6  - G       -- Soil heat flux, W/m-2
%              ::  7  - LAI     -- Leaf area index
%              ::  8  - soilpar -- Parameters related to Soil types
%              ::  9  - pftpar  -- Parameters related to PFTs
%              ::  10 - wa      -- Soil moisture (last step)
%              ::  11 - zgw     -- groundwater table depth, mm
%              ::  12 - snp     -- Snow package (old), mm day-1
% -------------------------------------------------------------------------
% Model output ::  1 - Et      -- Total Evapotranspiration, mm day-1
%              ::  2 - Tr      -- Plant Transpiration, mm day-1
%              ::  3 - Es      -- Soil Evaporation, mm day-1
%              ::  4 - Ei      -- Intercepted Evaporation, mm day-1
%              ::  4 - Esb     -- Snow sublimation, mm day-1
%              ::  5 - wa      -- Soil moisture (three layers)
%              ::  6 - srf     -- Surface runoff, mm day-1
%              ::  7 - zgw     -- groundwater table depth, mm
%              ::  8 - snp     -- Snow package (new), mm day-1
% -------------------------------------------------------------------------

% parameter section-----
alpha = optpara(1); % alpha = 1.26;
D50   = optpara(2);
D95   = optpara(3);
c     = -2.944/log(D95/D50);
% parameter section-----

% update pftpar
pftpar(2) = D50;
pftpar(3) = c;

% set the soil depth for three soil layers
zm = [50, 1450, 3500]; % mm

% potential Evaporation allocated to canopy and soil surface
[pEc, pEs] = potentialET(Rn, G, LAI, Ta, Pa, alpha);

% interception evaporation
[Ei, fwet, ~] = interception(LAI, Pe, pEc, pftpar);

% snow sublimation, snow melt
[snp, Esb, ~, Pnet] = snp_balance(Pe, Ta, Ei, snp);

% runoff
[srf, IWS, Vmax] = runoff_up(Pnet, zgw, zm, wa, soilpar);

% variables assciated with soil water balance
[wa, zgw, Tr, Es, uex] = sw_balance(IWS, pEc, pEs, Ta, Topt, s_VOD,...
    wa, soilpar, pftpar, fwet, zm, zgw); 

% total Evaporanspiration
Et = Tr + Es + Ei + Esb;
srf = srf + uex;

end
