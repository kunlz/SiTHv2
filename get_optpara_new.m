function [optpara] = get_optpara_new(LC)


if LC == 11
    Typenum = 13;
elseif LC == 22
    Typenum = 12;
elseif LC == 33
    Typenum = 10;
elseif LC == 40
    Typenum = 4;
elseif LC == 41
    Typenum = 1;
elseif LC == 42
    Typenum = 2;
elseif LC == 43
    Typenum = 3;
elseif LC == 44
    Typenum = 4;
elseif LC == 45
    Typenum = 5;
elseif LC == 55
    Typenum = 10;
elseif LC == 66
    Typenum = 10;
elseif LC == 77
    Typenum = 10;
elseif LC == 0
    Typenum = 10;
end


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

% % -------------- The optpar look up table -------------------------------01
% T =[0.98	771.6	2998.6; % 1     Evergreen Needleleaf Forest  ------ ENF
%     1.02	687.46	2816.1; % 2     Evergreen Broadleaf  Forest  ------ EBF
%     0.94	668.94	2936.1; % 3     Deciduous Needleleaf Forest  ------ DNF
%     0.96	1001.5	3188.4; % 4     Deciduous Broadleaf  Forest  ------ DBF
%     0.98	748.79	2866.9; % 5     Mixed                Forest  ------ MF
%     1.14	501.41	2163.2; % 6     Closed Shrubland             ------ CSH
%     1.08	753.08	2872.1; % 7     Open   Shrubland             ------ OSH
%     1.10	565.44	2753.4; % 8     Woody  Savannas              ------ WSA
%     1.06	759.28	2877.4; % 9     Savannas                     ------ SAV
%     1.26	598.88	2523.7; % 10    Grassland                    ------ GRA
%     1.30	896.84	3103.6; % 11    Cropland                     ------ CRO
%     1.12	752.93	2809.6];% 12    Wetland                      ------ WET

% % -------------- The optpar look up table -------------------------------02
% T =[0.98	771.6	2998.6; % 1     Evergreen Needleleaf Forest  ------ ENF
%     0.98	687.46	2816.1; % 2     Evergreen Broadleaf  Forest  ------ EBF
%     0.94	668.94	2936.1; % 3     Deciduous Needleleaf Forest  ------ DNF
%     0.96	1001.5	3188.4; % 4     Deciduous Broadleaf  Forest  ------ DBF
%     0.98	748.79	2866.9; % 5     Mixed                Forest  ------ MF
%     1.14	501.41	2163.2; % 6     Closed Shrubland             ------ CSH
%     1.12	753.08	2872.1; % 7     Open   Shrubland             ------ OSH
%     1.08	565.44	2753.4; % 8     Woody  Savannas              ------ WSA
%     1.07	759.28	2877.4; % 9     Savannas                     ------ SAV
%     1.26	598.88	2523.7; % 10    Grassland                    ------ GRA
%     1.32	896.84	3103.6; % 11    Cropland                     ------ CRO
%     1.22	752.93	2809.6];% 12    Wetland                      ------ WET

% -------------- The optpar look up table -------------------------------03
T =[1.26	1771.6	3998.6; % 1     Evergreen Needleleaf Forest  ------ ENF
    1.26	2187.4	4316.1; % 2     Evergreen Broadleaf  Forest  ------ EBF
    1.26	1668.9	3936.1; % 3     Deciduous Needleleaf Forest  ------ DNF
    1.26	1401.5	3888.4; % 4     Deciduous Broadleaf  Forest  ------ DBF
    1.26	1748.7	3866.9; % 5     Mixed                Forest  ------ MF
    1.26	701.41	2663.2; % 6     Closed Shrubland             ------ CSH
    1.26	753.08	2872.1; % 7     Open   Shrubland             ------ OSH
    1.26	765.44	2253.4; % 8     Woody  Savannas              ------ WSA
    1.26	559.28	1877.4; % 9     Savannas                     ------ SAV
    1.26	598.88	1823.7; % 10    Grassland                    ------ GRA
    1.26	896.84	2203.6; % 11    Cropland                     ------ CRO
    1.26	752.93	2809.6];% 12    Wetland                      ------ WET

% -------------- The optpar look up table ---------- globmap   ------------
% T =[0.98	771.6	2998.6; % 1     Evergreen Needleleaf Forest  ------ ENF
%     1.02	687.46	2816.1; % 2     Evergreen Broadleaf  Forest  ------ EBF
%     0.94	668.94	2936.1; % 3     Deciduous Needleleaf Forest  ------ DNF
%     0.96	1001.5	3188.4; % 4     Deciduous Broadleaf  Forest  ------ DBF
%     0.98	748.79	2866.9; % 5     Mixed                Forest  ------ MF
%     1.14	501.41	2163.2; % 6     Closed Shrubland             ------ CSH
%     1.14	753.08	2872.1; % 7     Open   Shrubland             ------ OSH
%     1.10	565.44	2753.4; % 8     Woody  Savannas              ------ WSA
%     1.12	759.28	2877.4; % 9     Savannas                     ------ SAV
%     1.30	598.88	2523.7; % 10    Grassland                    ------ GRA
%     1.32	896.84	3103.6; % 11    Cropland                     ------ CRO
%     1.16	752.93	2809.6];% 12    Wetland                      ------ WET



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