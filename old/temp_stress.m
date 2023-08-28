% ---------------------- %
% Temperature Constrains %
% ---------------------- %
function [St] = temp_stress(Topt,Ta)
% ------ function input -------
% Topt  : Optimum temperature for plant growning
% Ta    : Air temperature
% ------ function output -------
% St    : Temperature stress
% ------
St = exp(-((Ta-Topt)./Topt).^2);
end