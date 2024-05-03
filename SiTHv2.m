% ----------------------------------------------------------------------- %
%                                                                         %
%                            ==============                               %
%                             SiTHv2 model                                %
%                            ==============                               %
%                                                                         %
% Written by : Kun Zhang & Gaofeng Zhu                                    %
%              HKU, LZU                                                   %
% Original version : 5/9/2018                                             %
% Last modified : 17/1/2023                                               %
% Questions to zhangkun322@foxmail.com                                    %
% ----------------------------------------------------------------------- %

function [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, Pnet, IWS, Vmax] = SiTHv2(Rn, Ta, Tas, Topt, ...
    Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp)
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

% set the soil depth for three soil layers
zm = [50, 1450, 3500]; % mm

% potential Evaporation allocated to canopy and soil surface
[pEc, pEs] = potentialET(Rn, G, LAI, Ta, Pa);

% interception evaporation
[Ei, fwet, ~] = interception(LAI, Pe, pEc, pftpar);

% snow sublimation, snow melt
new_Pe = max(Pe - Ei, 0);
[snp, Esb, ~, Pnet] = snp_balance(new_Pe, Ta, Tas, snp, pEs);

% runoff
[srf, IWS, Vmax] = runoff_up(Pnet, zgw, zm, wa, soilpar);

% variables assciated with soil water balance
new_pEs = max(pEs - Esb, 0);
[wa, zgw, Tr, Es, uex] = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, ...
    wa, soilpar, pftpar, fwet, zm, zgw);

% total Evaporanspiration
Et = Tr + Es + Ei + Esb;
srf = srf + uex;

end
