% -------------------------- %
%  Interception evaporation  %
% -------------------------- %
function [Ei,wetT] = interception(pEc,LAI,Rain,PFTpar3)
% ------- function input -------
% pEc    : potnetial Evaporation on canopy
% Rain   :
% LAI    : Leaf area index
% PFTpar : PFT parameters
% ------- function output ------
% Ei     : Interception envaporation, mm/day
% wet    : wetness fraction
% -------
% Reference:
% Choudhury BJ, Digirolamo NE, 1998, A biologisical process-based estimate
% of global land surface evaporation using satellite and ancillary data I.
% Model description and comparison with observation. Journal of Hydrology.
%
% Interception (I, mm) was calculated using Horton's model:
%   I       =   min(P,aP+b)
% Where P   =   precipitation in mm
%  a,b      =   parameters with values derived from published records
% -------------------------------------------------------------------------

% x, parameter values, used for calculating daily rainfall interception
inc = PFTpar3;
SI = min(Rain,inc.*LAI.*Rain);
wetT = min(x.*SI./pEc,1);

if pEc < 1e-3

    Ei = 0;
else
    Ei = pEc.*wetT; 
end

end