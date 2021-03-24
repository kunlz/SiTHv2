function soilpar = SoilPara(SC)
% ----------Input----------
% SC      :: soilcode raster 
% ----------Output----------
% soilpar ::
%     1  ks 
%     2  theta_sat 
%     3  soil layer 1 water holding capacity (%) |
%     4  soil layer 2 water holding capacity (%) |
%     5  soil layer 2 water holding capacity (%) |
%     6  theta_wp 
%     7  theta_c
%     8  Wc
%     9  theta_FC

% ----------- 3-D rasters that contians the seven soil characteristic parameters ---------------
%        1           2         3         4       5         6        7
%        Ksat   theta_sat  theta_sat	 b    theta_FC  theta_c	 theta_wp  Soil Type
Store = [0.0236	   -47       0.373	   3.39	   0.151	 0.109	  0.035    % 1 Sand
         0.0166	   -64       0.386	   3.86	   0.189	 0.142	  0.052    % 2 Loamy sand
         0.0071	   -132	     0.416	   4.50	   0.265	 0.208	  0.087    % 3 Sandy loam
         0.0042	   -207	     0.435	   5.77	   0.331	 0.274	  0.139    % 4 Loam
         0.0042    -207      0.435     5.77    0.331     0.274    0.139    % 5 Slit
         0.0017	   -454	     0.468	   4.98	   0.399	 0.320	  0.146    % 6 Silty loam
         0.0071	   -132	     0.416	   7.20	   0.314	 0.270	  0.157    % 7 Sandy clay loam
         0.0028	   -289	     0.449	   8.32	   0.387	 0.339	  0.212    % 8 Clay loam
         0.0013	   -561	     0.476	   8.32	   0.444	 0.389	  0.243    % 9 Silty clay loam
         0.0058	   -158	     0.423	   9.59	   0.349	 0.312	  0.207    % 10 Sandy clay
         0.0011	   -633	     0.481	   10.4	   0.460	 0.414	  0.284    % 11 Silty clay
         0.0020	   -391	     0.461	   12.1	   0.427	 0.390	  0.282];  % 12 clay
% -----------------------------------------------------------------------------------------------
[m,n] = size(SC);
% tem = unique(SC,'sorted');
SPLUT = zeros(m,n,7);

for k = 1 : 7
    A = zeros(m,n);
    for i = 1 : 12
       A(SC == i) = Store(i,k);  
    end
    SPLUT(:,:,k) = A;
end
% Define ks value (texture-dependent)
soilpar(:,:,1) = SPLUT(:,:,1).*24.*3600; % mm day-1

% Define theta_sat value (texture-dependent)
soilpar(:,:,2) = SPLUT(:,:,2);

% Define WHCs of the 3 soil layers - Soil is uniformly distributed 
whc = SPLUT(:,:,5) - SPLUT(:,:,7);
soilpar(:,:,3) = whc;
soilpar(:,:,4) = whc;
soilpar(:,:,5) = whc;

% Define the wilting point
soilpar(:,:,6) = SPLUT(:,:,7);

% Define the critical values
soilpar(:,:,7) = SPLUT(:,:,6);

% wc
soilpar(:,:,8) = SPLUT(:,:,6)-SPLUT(:,:,7);

% theta_fc
soilpar(:,:,9) = SPLUT(:,:,5);
end