% ------------------------ %
%     Snowpack balance     %
% ------------------------ %
function [snowpack, Esb, snowmelt, Pnet] = snp_balance(preci, Ta, Tas, snowpack, pEs)
% snowpack :: available snow storage
% snowmelt :: snow melt
% Esb      :: snow sublimation

% Esnow_emp = 0.84.*(0.864 .* (7.093 .* Ta + 28.26)) ./ (Ta.^2 - 3.593 .* Ta + 5.175);
Esnow = pEs; % Simple equivalent (Needs further development)

% only snowfall occurs at Ta below zero
if Ta <= 0
    
    % Add new snowfall, Ta<=0
    newsnow = preci;
    snowpack = snowpack + newsnow;
    
    % snon melt
    snowmelt = 0;
    
    % real snow sublimation
    Esb = min(snowpack, Esnow);
    Esb = max(Esb, 0); % >0
    
    % net Precipitation into soil surface
    Pnet = 0;

    % new snowpack
    snowpack = snowpack - Esb;
    
else
    
    % real snow sublimation
    Esb = min(snowpack, Esnow);
    Esb = max(Esb, 0);

    snowpack = snowpack - Esb;

    % snow melt, Ta>0
    snowmelt_x = (1.5 + 0.007 * preci) * Tas; % Tas, accumulated Ta > 0
    snowmelt = min(snowpack, snowmelt_x);
    snowpack = snowpack - snowmelt;
    
    % net water into soil surface
    Pnet = max(0, preci + snowmelt);
    
end

end
