% --------------- %
%      Runoff     %
% ----------------%
function [srf, d] = runoff(Pnet,z_gwt,soildeep,swc,soilpar)
% --------- function input -------
% SoilDeep : Soil depth of different layers
% wa       : The antecedent soil water content expressed
%            as a function of the WHC in that layer
% soilpar  : Soil parameters according to Soil type
% Prnet    : Net precipitation = P-I+Snowmelt
% --------- function output ------
% Srf      : Surface Runoff
% d        : Thickness of the unsaturated soil in the ith layer
% ---------
% Reference:
% Choudhury BJ, Digirolamo NE, 1998
% SCS (1985). National Engineering Handbook. Section 4: Hydrology.
% Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.
% SCS (1986). Urban Hydrology for Small Watersheds, Technical Release No. 55.
% Washington, DC: Soil Conservation Service, U.S. Department of Agriculture.
% -------------------------------------------------------------------------

% Soil layers numbers 
lynum = length(soildeep);

% water holding capacity 
whc = soilpar(3:5); 

zm1 = soildeep(1);
zm2 = soildeep(2);
zm3 = soildeep(3);

% calculate the thickness of unsaturated soil in ith layer, (mm)
% zgw = groundwater table depth 
if z_gwt <= 0
    d1 = 0;
    d2 = 0;
    d3 = 0;
    
elseif z_gwt>0 && z_gwt<=zm1 
    d1 = z_gwt; 
    d2 = 0;
    d3 = 0;
    
elseif z_gwt>zm1 && z_gwt<=zm2
    d1 = zm1;
    d2 = z_gwt-zm1; 
    d3 = 0;
    
elseif z_gwt>zm2 && z_gwt<=zm3
    d1 = zm1;
    d2 = zm2;
    d3 = z_gwt-zm2-zm1;
    
else
    d1 = zm1;
    d2 = zm2;
    d3 = zm3;
    
end

d = [d1,d2,d3]; 

% calculate the overall soil water retention capacity
Vmax = 0;
for i = 1 : lynum 
    Vmax = Vmax+(1-swc(i)).*d(i).*whc(i);
end

% calculate the surface runoff, mm day-1
if Pnet > 0.2*Vmax
    srf = (Pnet-Ia1).^2./(Pnet+0.8*Vmax);
else
    srf = 0;
end

end