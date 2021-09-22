function [optpara] = get_optpara(Typenum)

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

% -------------- The optpar look up table ---------------------------------
T =[0.98	771.6	2998.6; % 1     Evergreen Needleleaf Forest  ------ ENF
    1.02	687.46	2816.1; % 2     Evergreen Broadleaf  Forest  ------ EBF
    0.94	668.94	2936.1; % 3     Deciduous Needleleaf Forest  ------ DNF
    0.96	1001.5	3188.4; % 4     Deciduous Broadleaf  Forest  ------ DBF
    0.98	748.79	2866.9; % 5     Mixed                Forest  ------ MF
    1.14	501.41	2163.2; % 6     Closed Shrubland             ------ CSH
    1.08	753.08	2872.1; % 7     Open   Shrubland             ------ OSH
    1.10	565.44	2753.4; % 8     Woody  Savannas              ------ WSA
    1.06	759.28	2877.4; % 9     Savannas                     ------ SAV
    1.26	598.88	2523.7; % 10    Grassland                    ------ GRA
    1.30	896.84	3103.6; % 11    Cropland                     ------ CRO
    1.12	752.93	2809.6];% 12    Wetland                      ------ WET


% typeclass = {'ENF','EBF','DNF','DBF','MF','CSH','OSH','WSA','SAV','GRA','CRO','WET'};
% plantType = 'MF';
% [~,b] = ismember(plantType, typeclass);

if Typenum == 0 || Typenum == 13 || Typenum == 15 
    optpara = T(10, :); % set to 10
    
elseif Typenum == 11
    optpara = T(12, :); % set to Wetland
    
elseif Typenum == 12 || Typenum == 14 || Typenum == 16
    optpara = T(11, :); % set to Cropland
    
else
    optpara = T(Typenum, :); % set to Wetland
    
end

end