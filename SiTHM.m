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


function [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp] = SiTHM(Rn, Ta, Topt,...
    Pe, Pa, G, LAI, soilpar, pftpar, wa, zgw, snp)
% Main function
% -------------------------------------------------------------------------
% Model inout  ::  1 - Rn      -- Net Radiation, W/m-2
%              ::  2 - Ta      -- Near surface air temperature, C
%              ::  2 - Topt    -- Optimal growth temperature for plant, C
%              ::  3 - Pe      -- Precipetition, mm day-1
%              ::  4 - Pa      -- Near surface air pressure, kPa
%              ::  5 - G       -- Soil heat flux, W/m-2
%              ::  6 - LAI     -- Leaf area index
%              ::  7 - soilpar -- Parameters related to Soil types
%              ::  8 - pftpar  -- Parameters related to PFTs
%              ::  9 - wa      -- Soil moisture (last step)
%              ::  10- zgw     -- groundwater table depth, mm
%              ::  11- snp    -- Snow package (old), mm day-1
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
zm = [50, 550, 2400]; % mm

% potential Evaporation allocated to canopy and soil surface
[pEc, pEs] = potentialET(Rn, G, LAI, Ta, Pa);

% interception evaporation
[Ei, fwet] = interception(LAI, Pe, pEc, pftpar);

% snow sublimation, snow melt
[snp, Esb, snm] = snp_balance(Pe, Ta, snp);

% net Precipitation into soil surface
Pnet = max(0, Pe + snm - Esb - Ei); 

% runoff
[srf, IWS] = runoff(Pnet, zgw, zm, wa, soilpar);

% variables assciated with soil water balance
[wa, zgw, Tr, Es, uex] = sw_balance(IWS, pEc, pEs, Ta, Topt, wa,...
    soilpar, pftpar, fwet, zm, zgw);

% total Evaporanspiration
Et = Tr + Es + Ei + Esb;
srf = srf + uex;

end
