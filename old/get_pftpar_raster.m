function [pftpar] = get_pftpar_raster(Typenum)

% ----------Input----------
% Plant Type -- according to MCD12, IGBP
% 0     'water'
% 1     'evergreen needleleaf forest'
% 2     'evergreen broadleaf forest'
% 3     'deciduous needleleaf forest'
% 4     'deciduous broadleaf forest'
% 5     'mixed forests'
% 6     'closed shrubland'
% 7     'open shrublands'
% 8     'woody savannas'
% 9     'savannas'
% 10	'grasslands'
% 11	'permanent wetlands'
% 12	'croplands'
% 13	'urban and built-up'
% 14	'cropland/natural vegetation mosaic'
% 15	'snow and ice'
% 16	'barren or sparsely vegetated'

% ----------Output---------
% pftpar  ::
%  1  beta -- the coefficient for calculation of interceptions
%  2  D50  -- the depth above which 50% of the root mas is located, mm
%  3  c    -- the shape parameters of logistic dose-response root distribution model
%  4  Zr   -- root depth (m)

% -------------- The PFTpar look up table ---------------------------------
%       beta   D50      c     Zr    CH   %       PFT
PLUT = [0.06   120    -1.87   2.3   25;     % 1     Evergreen Needleleaf Forest   ------ ENF
        0.02   200    -1.67   7.3   20;     % 2     Evergreen Broadleaf  Forest   ------ EBF
        0.06   120    -1.87   2.0   20;     % 3     Deciduous Needleleaf Forest   ------ DNF
        0.06   120    -1.87   2.0   20;     % 4     Deciduous Broadleaf  Forest   ------ DBF
        0.04   180    -1.87   2.8   15;     % 5     Mixed                Forest   ------ MF
        0.01   190    -1.34   4.8   2;      % 6     Closed Shrubland              ------ CSH
        0.01   280    -1.91   5.2   2;      % 7     Open   Shrubland              ------ OSH
        0.01   280    -1.80   5.0   2;      % 8     Woody  Savannas               ------ WSA
        0.01   230    -1.63   8.0   0.5;    % 9     Savannas                      ------ SAV
        0.01   130    -1.92   2.6   0.5;    % 10    Grassland                     ------ GRA
        0.01   180    -2.01   2.1   2;      % 11    Cropland                      ------ CRO
        0.01   120    -1.90   2.3   2];     % 12    Wetland                       ------ WET
% -------------- The PFTpar look up table ---------------------------------


% typeclass = {'ENF','EBF','DNF','DBF','MF','CSH','OSH','WSA','SAV','GRA','CRO','WET'};
% plantType = 'MF';
% [~,b] = ismember(plantType, typeclass);

if Typenum == 0 || Typenum == 13 || Typenum == 15 || Typenum == 16
    pftpar = PLUT(10, :); % set to 10
    
elseif Typenum == 11
    pftpar = PLUT(12, :); % set to Wetland
    
elseif Typenum == 12 || Typenum == 14
    pftpar = PLUT(11, :); % set to Cropland
    
else
    pftpar = PLUT(Typenum, :); % set to Wetland
    
end


end