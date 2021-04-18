function [pftpar] = get_pftpar(plantType)

% ----------Input----------
% Plant Type 
% ----------Output---------
% pftpar  ::
%  1  beta -- the coefficient for calculation of interceptions
%  2  D50  -- the depth above which 50% of the root mas is located, mm
%  3  c    -- the shape parameters of logistic dose-response root distribution model
%  4  Zr   -- root depth (m)

% -------------- The PFTpar look up table ---------------------------------
%       beta   D50      c     Zr       %       PFT
PLUT = [0.06   120	  -1.87   2.3;     % 1     Evergreen Needleleaf Forest   ------ ENF
        0.02   200    -1.67   7.3;     % 2     Evergreen Broadleaf  Forest   ------ EBF
        0.06   120	  -1.87   2.0;     % 3     Deciduous Needleleaf Forest   ------ DNF
        0.06   120	  -1.87   2.0;     % 4     Deciduous Broadleaf  Forest   ------ DBF
        0.04   180	  -1.87   2.8;     % 5     Mixed                Forest   ------ MF
        0.01   190    -1.34   4.8;     % 6     Closed Shrubland              ------ CSH
        0.01   280    -1.91   5.2;     % 7     Open   Shrubland              ------ OSH
        0.01   280    -1.80   15.0;    % 8     Woody  Savannas               ------ WSA
        0.01   230    -1.63   8.0;     % 9     Savannas                      ------ SAV
        0.01   130    -1.92   2.6;     % 10    Grassland                     ------ GRA
        0.01   180    -2.01   2.1;     % 11    Cropland                      ------ CRO
        0.01   120    -1.90   2.3];    % 12    Wetland                       ------ WET
% -------------- The PFTpar look up table ---------------------------------
pftpar = PLUT(plantType, :);

end


