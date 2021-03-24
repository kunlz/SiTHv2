function [Rnet] = cal_Rn(SWD, Ta, albedo)

% RH = SH2RH(SH,Ta,Pa);
% RH(RH<0) = 0;
% RH(RH>100) = 100;

Rnet = ISW2NET(SWD,Ta-273.16,albedo);

end

% function [RH] = SH2RH(SH,T,Pa)
% % note:
% % this function designed for calculate Net Radiation
% % input: SH   : specific humidity, 1
% %        T    : Air temperature, K
% %        Pa   : air pressure, Pa
% RH = 0.263.*Pa.*SH.*(exp(17.67.*(T-273.16)./(T-29.65))).^(-1);
% end

function [Rnet] = ISW2NET(ISW,T,alfa)
% note:
% this function designed for calculate Net Radiation
% input: SW_in: incoming shortwave radiation
%        T    : Air temperature in C
%        alfa : MODIS albedo
% alfa(alfa>1000) = 400;
alfa = 0.0001.*alfa;
alfa(alfa<0) = 0;
sigma = 5.678e-8;
epsilong_s = 0.97;
epsilong_a = 1 - 0.26.*exp(-7.77 .* 1e-4 .* T.^2);
Rnet = (1-alfa).*ISW + (epsilong_a-epsilong_s).*sigma.*(273.15+T).^4;
% Rnet(Rnet<0) = 0;
end