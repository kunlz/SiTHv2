% ----------------------------------------------------------------------- %                           
%                                                                         %
%                            ==============                               %
%                              STHM model                                 %
%                            ==============                               %
%                                                                         %                          
% Written by : Kun Zhang, Gaofeng Zhu;  Lanzhou University, China         %
% Date       : 5/9/2018                                                   %
% Questions to zhangkun322@foxmail.com                                    %
% ----------------------------------------------------------------------- %
% Last modified : 16/3/2021, vector version

function [Et,Tr,Es,Ei,Esnow,wa,Srf,SnpN] = STHM(Rn,Ta,Topt,Rain,Pa,G,LAI,PFTpar,soilpar,wo,SnpO)
% Main function
% -------------------------------------------------------------------------
% Model inout  ::  1 - Rn      -- Net Radiation, W/m-2
%              ::  2 - Ta      -- Near surface air temperature, C
%              ::  3 - Rain    -- Precipetition, mm/day
%              ::  4 - Pa      -- Near surface air pressure, kPa
%              ::  5 - G       -- Soil heat flux, W/m-2
%              ::  6 - LAI     -- Leaf area index
%              ::  7 - PFTpar  -- Parameters related to PFTs
%              ::  8 - Soilpar -- Parameters related to Soil types
%              ::  9 - wo      -- Soil moisture (last step)
%              :: 10 - SnpN    -- Snow package (old),  mm/day
% -------------------------------------------------------------------------
% Model output ::  1 - Et      -- Total Evapotranspiration, mm/day
%              ::  2 - Tr      -- Plant Transpiration, mm/day
%              ::  3 - Es      -- Soil Evaporation, mm/day
%              ::  4 - Ei      -- Intercepted Evaporation, mm/day
%              ::  4 - Esnow   -- Snow sublimation, mm/day
%              ::  5 - wa      -- Soil moisture (three layers), %
%              ::  6 - Srf     -- Surface runoff,  mm/day
%              ::  7 - SnpN    -- Snow package (new),  mm/day
% -------------------------------------------------------------------------

% Set the soil depth for three soil layers
soildeep = [50,500,2000]; % mm

% Potential Evaporation allocated to canopy and soil surface
[pEc,pEs] = potentialET(Rn,G,LAI,Ta,Pa);

% Interception evaporation
[Ei,wet] = interception(LAI,Rain,pEc,PFTpar);

% Snow sublimation
[SnpN,Esnow] = snowpack_balance(Rain,Ta,SnpO);

% Calculate the net Precipitation into soil surface
Prnet = Rain-Esnow-Ei;
Prnet(Prnet<0) = 0;

% Calculate the Runoff
[Srf] = Runoff(soildeep,wo,soilpar,Prnet);

% Partition the potential tranpiration to three layers
[pTr_ly] = pTrPartition(pEc,wo,PFTpar,wet);




% Calculate the Infiltration into soil surface
IWS = Prnet-Srf; 

% Calculate Plant transpiration, Soil evaporation, and Soil moisture status
[Tr,Es,wa,~,~] = SWB(Topt,Ta,IWS,pTr_ly,pEs,wo,wet,soilpar,soildeep);

% Calculate the total Evaporanspiration
Tr(Tr<0) = 0;
Es(Es<0) = 0;
Ei(Ei<0) = 0;
Esnow(Esnow<0) = 0;
Et = Tr+Es+Ei+Esnow;

end

