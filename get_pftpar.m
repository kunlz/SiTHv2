function [pftpar] = get_pftpar(LC)

% ---------------------------------------------------------------- %
% Standard numbering rule for different PFTs based on MCD12, IGBP- %
% ---------------------------------------------------------------- %
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
%  pftpar  ::
%  1  beta -- the coefficient for calculation of interceptions
%  2  D50  -- the depth above which 50% of the root mas is located, mm
%  3  c    -- the shape parameters of logistic dose-response root distribution model
%  4  Zr   -- root depth (m)

% -------------- The PFTpar look up table ---------------------------------
%       beta   D50    D95    CH      %       PFT
PLUT = [0.06   1771   3998   20;     % 1     Evergreen Needleleaf Forest   ------ ENF
        0.02   2187   4316   18;     % 2     Evergreen Broadleaf  Forest   ------ EBF
        0.06   1668   3936   18;     % 3     Deciduous Needleleaf Forest   ------ DNF
        0.06   1401   3888   18;     % 4     Deciduous Broadleaf  Forest   ------ DBF
        0.04   1748   3866   15;     % 5     Mixed                Forest   ------ MF
        0.01   701    2663   1.5;    % 6     Closed Shrubland              ------ CSH
        0.01   753    2872   1.5;    % 7     Open   Shrubland              ------ OSH
        0.01   765    2253   1.2;    % 8     Woody  Savannas               ------ WSA
        0.01   559    1877   1.0;    % 9     Savannas                      ------ SAV
        0.01   598    1823   1.0;    % 10    Grassland                     ------ GRA
        0.01   896    2203   1.5;    % 11    Cropland                      ------ CRO
        0.01   752    2809   1.5];   % 12    Wetland                       ------ WET
% -------------- The PFTpar look up table ---------------------------------

% -------------------------------------------------------------------------
% When use a longterm landcover dataset: HILDA +
% Adjust the PFT code to match the IGBP-LC

% 1) Land Use/Cover categories (states)
% 11 Urban
% 22 Cropland
% 33 Pasture
% 40: Forest (Unknown/Other)
% 41: Forest (Evergreen, needle leaf)
% 42: Forest ( Evergreen, broad leaf)
% 43: Forest (Deciduous, needle leaf)
% 44: Forest (Deciduous, broad leaf)
% 45: Forest (Mixed)
% 55 Grass/shrubland
% 66 Other land
% 77 Water

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
    Typenum = 11;
elseif LC == 0
    Typenum = 10;
end
% -------------------------------------------------------------------------

% IGBP <-----> PFTpar table
if Typenum == 0 || Typenum == 13 || Typenum == 15 || Typenum == 16
    pftpar = PLUT(10, :); % set to Grassland
elseif Typenum == 11
    pftpar = PLUT(12, :); % set to Wetland 
elseif Typenum == 12 || Typenum == 14
    pftpar = PLUT(11, :); % set to Cropland
else
    pftpar = PLUT(Typenum, :);
end

end