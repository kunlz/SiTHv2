% ---------------------------------- %
%  Potential transpiration partition %
% -----------------------------------%
function [pTr_ly] = pTr_partition(pEc,wo,PFTpar,wet)
% ------- function input -------
% pEc    :  Potnetial Evaporation on canopy
% w      :  Initialized values for soil moisture
% PFTpar :  PFT parameters
% ------- function output ------
% pTr_ly :  seperate potnetial Transpiration
% -------

% Get all avaliable water contents through root distribution
wr = PFTpar(1).*wo(1)+PFTpar(2).*wo(2)+(1-PFTpar(1)-PFTpar(2)).*wo(3);

% Root distribution adjusted by soil water content
beta1 = PFTpar(1).*wo(1)./wr;
beta2 = PFTpar(2).*wo(2)./wr;
beta3 = (1-PFTpar(1)-PFTpar(2)).*wo(3)./wr;

beta1(wr==0) = 0;
beta2(wr==0) = 0;
beta3(wr==0) = 0;

% Get the potentail Transpiration for each layer
pTr_ly(:,:,1) = beta1.*pEc.*(1-wet);
pTr_ly(:,:,2) = beta2.*pEc.*(1-wet);
pTr_ly(:,:,3) = beta3.*pEc.*(1-wet);

end