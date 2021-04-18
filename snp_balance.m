% ------------------------ %
%     Snowpack balance     %
% ------------------------ %
function [snowpack, Esb, snowmelt] = snp_balance(rain, Ta, snowpack)
% snowpack :: available snow storage
% snowmelt :: snow melt
% Esb      :: snow sublimation

% we set only snowfall occurs at Ta below zero
if Ta <= 0
    
    % Add new snowfall, Ta<=0
    newsnow = rain;
    snowpack = snowpack + newsnow;
    snowmelt = 0;
    
else
    
    % snow melt, Ta>0
    snowmelt_x = (1.5 + 0.007 * rain) * Ta;
    snowmelt = min(snowpack, snowmelt_x);
    snowpack = snowpack - snowmelt;
    
end

% snow sublimation, set as an emperical factor
Esnow_emp = (0.0864 * (-7.093 * Ta + 28.26)) / (Ta^2 - 3.593 * Ta + 5.175);
Esb = min(snowpack, Esnow_emp);

% new snowpack
snowpack = snowpack - Esb;

end
