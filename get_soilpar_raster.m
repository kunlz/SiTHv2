function [soilpar] = get_soilpar_raster(SC)
% ----------Input----------
% Soil Type idnex
% ----------Output----------
% soilpar ::
%     1  ks         -- hydraulic conductivity, m day-1
%     2  pha_sat    -- water potential at saturation condition
%     3  theta_sat  -- water content at saturation condition
%     4  b          -- an empirical parameter
%     5  theta_fc   -- soil moisture at field capacity
%     6  theta_c    -- soil moisture at critical value
%     7  theta_wp   -- soil moisture at wilting point
% --------------------

% ----------- 3-D rasters that contians the seven soil characteristic parameters ---------------
%        1            2            3        4        5          6       7
%        Ksat      pha_sat	   theta_sat    b     theta_fc  theta_c	 theta_wp   Soil Type
store = [1.3841	   -0.0232       0.373	   3.39	    0.151	  0.109	  0.035;    % 1 Sand
         0.8229	   -0.0175       0.386	   3.86	    0.189	  0.142	  0.052;    % 2 Loamy sand
         0.5353	   -0.0316	     0.416	   4.50	    0.265	  0.208	  0.087;    % 3 Sandy loam
         0.4086	   -0.0655	     0.435	   5.77	    0.331	  0.274	  0.139;    % 4 Loam
         0.4331	   -0.0562	     0.455	   5.32	    0.401	  0.317	  0.189;    % 5 Slit
         0.4427	   -0.0471	     0.468	   4.98	    0.399	  0.320	  0.146;    % 6 Silty loam
         0.4991	   -0.0310	     0.416	   7.20	    0.314	  0.270	  0.157;    % 7 Sandy clay loam
         0.3552	   -0.0599	     0.449	   8.32	    0.387	  0.339	  0.212;    % 8 Clay loam
         0.3848	   -0.0414	     0.476	   8.32	    0.444	  0.389	  0.243;    % 9 Silty clay loam
         0.6157	   -0.0269	     0.423	   9.59	    0.349	  0.312	  0.207;    % 10 Sandy clay
         0.2967	   -0.0453	     0.481	   10.4	    0.460	  0.414	  0.284;    % 11 Silty clay
         0.2580	   -0.0531	     0.461	   12.1	    0.427	  0.390	  0.282];   % 12 clay
% -----------------------------------------------------------------------------------------------
store(:,1) = 1000 .* store(:,1);

% soilty = {'Sand','Loamy_Sand','Sandy_Loam','Loam','Slit','Silty_Loam',...
%           'Sandy_Clay_Loam','Clay_Loam','Silty_Clay_Loam','Sandy_Clay',...
%           'Silty_Clay','Clay'};

% [~,b] = ismember(soilType, soilty);      
if SC == 0
    SC = 12;
end
soilpar = store(SC, :);

end
