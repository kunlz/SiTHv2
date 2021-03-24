function [PFTpar] = PFTPara(LC)
% ----------Input----------
% LC      :: Landcover Type raster
% ----------Output----------
% PFTpar  ::
%  1  z1 -- the fraction of the plants roots which are in the top soil layer
%  2  z2 -- the fraction of the plants roots which are in the second and third layers
%  3  i  -- canopy interception coefficient

% -------------------------------- The PFTpar look up table ---------------------------------
%        z1     z2      i         %       PFT
PLUT = [0.00   0.60	  0.06        % 1     Evergreen Needleleaf Forest   ------ ENF
        0.00   0.70   0.02        % 2     Evergreen Broadleaf  Forest   ------ EBF
        0.00   0.90	  0.06        % 3     Deciduous Needleleaf Forest   ------ DNF
        0.00   0.70	  0.02        % 4     Deciduous Broadleaf  Forest   ------ DBF
        0.00   0.70	  0.02        % 5     Mixed                Forest   ------ MF
        0.00   0.80   0.01        % 6     Closed Shrubland              ------ CSH
        0.00   0.80   0.01        % 7     Open   Shrubland              ------ OSH
        0.00   0.80   0.01        % 8     Woody  Savannas               ------ WSA
        0.00   0.80   0.01        % 9     Savannas                      ------ SAV
        0.00   0.80   0.01        % 10    Grassland                     ------ GRA
        0.00   0.80   0.01        % 11    Cropland                      ------ CRO
        0.00   0.80   0.01];      % 12    Wetland                       ------ WET
% -------------------------------- The PFTpar look up table ---------------------------------
%%
tem = unique(LC,'sorted');
[m,n] = size(LC);
PFTpar = zeros(m,n,3);
for k = 1 : 3
    A = zeros(m,n);
    for i = 1 : length(tem)
        % Select the water & snow & unclassified & filled value
        if tem(i)==0 || tem(i)==15 || tem(i)==254 || tem(i)==255
            continue
            
            % Grass: grassland, urban and built-up, barren or sparsely vegetated
        elseif tem(i)==10 || tem(i)==13 ||tem(i)==16
            A(LC==tem(i)) = PLUT(10,k);
            
            % Wetland
        elseif tem(i)==11
            A(LC==tem(i)) = PLUT(12,k);
            
            % Crop & Cropland/Natural vegetation mosaic
        elseif tem(i)==12 || tem(i)==14
            A(LC==tem(i)) = PLUT(11,k);
            
            % the rest landcover type
        else
            A(LC==tem(i)) = PLUT(tem(i),k);
        end
    end
    PFTpar(:,:,k) = A;
end

end

