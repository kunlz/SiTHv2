% ------------------------ %
%     Snowpack balance     %
% ------------------------ %
function [snowpack, Esb, snowmelt, Pnet] = snp_balance(preci, Ta, Ei, snowpack)
% snowpack :: available snow storage
% snowmelt :: snow melt
% Esb      :: snow sublimation

% snow sublimation, set as an emperical factor
Esnow_emp = 0.24.*(0.0864 .* (-7.093 .* Ta + 28.26)) ./ (Ta.^2 - 3.593 .* Ta + 5.175);

% we set only snowfall occurs at Ta below zero
if Ta <= 0
    
    % Add new snowfall, Ta<=0
    newsnow = preci;
    snowpack = snowpack + newsnow;
    
    % snon melt
    snowmelt = 0;
    
    % real snow sublimation
    Esb = min(snowpack, Esnow_emp);
    Esb = max(Esb, 0);
    
    % net Precipitation into soil surface
    Pnet = 0;
    
else
    
    % snow melt, Ta>0
    snowmelt_x = (1.5 + 0.007 * preci) * Ta;
    snowmelt = min(snowpack, snowmelt_x);
    snowpack = snowpack - snowmelt;
    
    % real snow sublimation
    Esb = min(snowpack, Esnow_emp);
    Esb = max(Esb, 0);

    % net Precipitation into soil surface
    Pnet = max(0, preci + snowmelt - Ei);
    
end

% new snowpack
snowpack = snowpack - Esb;

end
