% ------------------------ %
%     Snowpack balance
% ------------------------ %
function [snowpack,Esnow] = snowpack_balance(Rain,Ta,snowpack)
% snowpack :: available snow storage
% dmelt    :: snow melt
% Esnow    :: snow sublimation

% we set only snowfall occurs at Ta below zero
if Ta <= 0 
    
    % Add new snowfall, Ta<=0   
    newsnow = Rain;
    snowpack = snowpack+newsnow;
else
    
    % snow melt, Ta>0
    snowmelt_x = (1.5+0.007.*Rain).*Ta;
    snowmelt = min(snowpack,snowmelt_x);
    snowpack = snowpack-snowmelt;
end

% snow sublimation, set as an emperical factor
Esnow_emp = (0.0864.*(-7.093.*Ta+28.26))./(Ta.^2-3.593.*Ta+5.175);
Esnow = min(snowpack,Esnow_emp);


% new snowpack
snowpack = snowpack-Esnow;

end

