
clc
clear 

% Set spatial resolution
res = 0.1; % degrees, global
r_m = 180 / res;
r_n = 360 / res;
% Output files folder :
Output_folder = 'Y:/SiTHv2_out_longterm/';

if exist(Output_folder, 'dir') == 0.
    mkdir(Output_folder);
end

% Forcing data folder :
Forcing_folder = 'Y:/ModelForcingData/01deg/'; 

% Sub-folders of Focing data
subfolder_Rn = 'Rn/RN_ERA5L/RN.ERA5L.GLASS.S100.A';
subfolder_Preci = 'Preci/ERA5L/P.ERA5L.scale001.A';
subfolder_Ta = 'Ta/Ta.MSWX.scale001.A';
subfolder_Pa = 'Pa/Pa.MSWX.scale001.A';
subfolder_LAI = 'LAI/GEOV2_LAI/THEIA_GEOV2_R01_AVHRR_LAI_A'; 
subfolder_LC = 'LC/hildap_vGLOB1.0f/hildap_vGLOB.A';
subfolder_VOD = 'Vegetation_Optical_Depth/VODCA/0.1deg/VODCA_Xband_A';

% Set georeference, for 0.25 degrees, global
latlim = [-90,90];
lonlim = [-180,180];
rasterSize = [720,1440];
RA = georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');

% Load soil type
Soilraster = load('inpara\Soilraster.mat');
Soilraster = Soilraster.Soilraster; 

% Load land mask
maskland = load('inpara\landmask01.mat');
maskland = uint8(maskland.mask2);

% Load the optimal temperature for plant growth
Topt = load('inpara\Topt.mat');
Topt = single(Topt.Topt_new);

% Parallel calculation
% parpool('local', 40);

% Main loops
for yr = 1982 : 2020

    clear waa zgww snpp
    disp(' ')
    disp(['------------------------- Start calculation for ' num2str(yr)])
    % spin-up year : 100 years
    if yr == 1981

        spinfg = 1; % need spin-up
        disp(' ')
        disp('start year ... spin-up ... set spinfg = 1')

        % Initialization
        waa = 0.25 .* ones(r_m, r_n, 3); % initial value for swc
        zgww = 5050 .* ones(r_m, r_n); % initial value for groundwater table
        snpp = zeros(r_m, r_n); % initial value for snowpack depth

        % load the updated variables
        % uptval = load('upt_vals_01_backup.mat');
        % uptval = uptval.X_upt;
    else
        spinfg = 0;
        disp(' ')
        disp('normal year ... set spinfg = 0')
        % load the updated variables
        uptval = load('upt_vals_01.mat');
        uptval = uptval.X_upt;

        waa = uptval(:, :, 1:3); 
        zgww = uptval(:, :, 4); 
        snpp = uptval(:, :, 5); 
    end

    % ----------------- %
    % Load Forcing Data %
    % ----------------- %
    % --------------------------------------------------------------------%
    % Net Radiation, W/m2
    disp(['Load Net Radiation ... For the year :: ' num2str(yr)])
    % EMO_Rn = matfile([Forcing_folder subfolder_Rn num2str(yr) '.mat']);
    EMO_Rn = load([Forcing_folder subfolder_Rn num2str(yr) '.mat']);
    EMO_Rn = EMO_Rn.RN; 

    % Air Temperature, C, 2m
    disp(['Load Air Temperature ... For the year :: ' num2str(yr)])
    % EMO_Ta = matfile([Forcing_folder subfolder_Ta num2str(yr) '.mat']);
    EMO_Ta = load([Forcing_folder subfolder_Ta num2str(yr) '.mat']);
    EMO_Ta = EMO_Ta.Ta; 

    % Precipitation, mm
    disp(['Load Precipitation ... For the year :: ' num2str(yr)])
    % EMO_Preci = matfile([Forcing_folder subfolder_Preci num2str(yr) '.mat']);
    EMO_Preci = load([Forcing_folder subfolder_Preci num2str(yr) '.mat']);
    EMO_Preci = EMO_Preci.P; 

    % Air Pressure, kPa
    disp(['Load Air Pressure ... For the year :: ' num2str(yr)])
    % EMO_Pa = matfile([Forcing_folder subfolder_Pa num2str(yr) '.mat']);
    EMO_Pa = load([Forcing_folder subfolder_Pa num2str(yr) '.mat']);
    EMO_Pa = EMO_Pa.Pa;

    % Satellite-based LAI
    disp(['Load Satellite-based LAI  ... For the year :: ' num2str(yr)])
    % EMO_LAI = matfile([Forcing_folder subfolder_LAI num2str(yr) '.mat']);
    EMO_LAI = load([Forcing_folder subfolder_LAI num2str(yr) '.mat']);
    EMO_LAItime = EMO_LAI.tt; 
    EMO_LAI = EMO_LAI.LAIx; 
    
    % Satellite-based VOD
    disp('Load Satellite-based VOD  ...')
    % EMO_VOD = matfile([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    EMO_VOD = load([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    EMO_VOD = EMO_VOD.VODCAy;

    % Satellite-based Landcover/PFTs
    disp(['Load Satellite-based Landcover  ... For the year :: ' num2str(yr)])
    LC_year = load([Forcing_folder subfolder_LC num2str(yr) '.mat']);
    LC_year = LC_year.LULC; 
    % --------------------------------------------------------------------%
    
    % Days of selected year
    days = yeardays(yr);


    % ------------------ %
    % Parallel Computing %
    % ------------------ %

    disp('Preallocate memory to each variables ... ')
    % 10 variables
    X_ET  = zeros(r_m, r_n, days,'double');
    X_Tr  = zeros(r_m, r_n, days,'double');
    X_Es  = zeros(r_m, r_n, days,'double');
    X_Ei  = zeros(r_m, r_n, days,'double');
    X_Esb = zeros(r_m, r_n, days,'double');
    X_SM1 = zeros(r_m, r_n, days,'double');
    X_SM2 = zeros(r_m, r_n, days,'double');
    X_SM3 = zeros(r_m, r_n, days,'double');
    X_RF  = zeros(r_m, r_n, days,'double');
    X_GW  = zeros(r_m, r_n, days,'double');
    
    disp('Start calculation ... ')
    X_upt = zeros(r_m, r_n, 5);

    ppm = ParforProgressbar(r_m, 'showWorkerProgress', false);
    parfor i = 1 : r_m
        
        % read each row, (latitude)
        Rnix = permute(EMO_Rn(i, :, :), [3, 2, 1]);
        Taix = permute(EMO_Ta(i, :, :), [3, 2, 1]);
        Precix = permute(EMO_Preci(i, :, :), [3, 2, 1]);
        Paix = permute(EMO_Pa(i, :, :), [3, 2, 1]);
        LAIix = permute(EMO_LAI(i, :, :), [3, 2, 1]);
        VODix = permute(EMO_VOD(i, :, :), [3, 2, 1]);

        % X_mat for different output variables
        X_vals = zeros(days, r_n, 10); % total
        X_upti = zeros(1, r_n, 5); % intermediate

        for j = 1 : r_n
            
            % check landmask
            if maskland(i, j) == 0
                continue
            end

            % update variables
            wa = reshape(waa(i, j, :), [1, 3]);
            zgw = zgww(i, j);
            snp = snpp(i, j);

            % Meteo forcing for each pixel, need rescale
            Rni = 0.01.*double(Rnix(:, j));
            Tai = 0.01.*double(Taix(:, j));
            Precii = 0.01.*double(Precix(:, j));            
            Pai = 0.01.*double(Paix(:, j));

            % Cal Tas
            Tasi = Tai;
            Tasi(Tasi < 0) = 0;
            Tasi = cumsum(Tasi);

            % Satellite LAI for each pixel
            LAIi = 0.01.*double(LAIix(:, j));
            xo = day(EMO_LAItime,"dayofyear");
            xi = 1:1:days;
            LAIii = interp1(xo', LAIi, xi', 'pchip', 'extrap');
            LAIii(LAIii < 0) = 0;

            % Cal G_soil, % Choudhury et al., 1987
            Gi = 0.4 .* Rni .* exp(-0.5 .* LAIii);

            % Cal VOD-stress
            VODi = 0.001.*double(VODix(1:days, j));
            % VODi = smooth(VODi, 7, 'moving');
            VODi(VODi < 0) = 0;
            s_VODi = (VODi ./ max(VODi)).^0.5;

            % Topt
            Top = Topt(i, j); 
            if isnan(Top)
                Top = 25; % Topt = 25 when Topt is NaN;
            end

            % Parameter-set for plant and soil
            % 1- PlantType
            PFTi = LC_year(i, j);
            pftpar = get_pftpar(PFTi);

            % 2- SoilType
            SC = Soilraster(i, j);
            soilpar = get_soilpar_raster(SC);

            % check wa
            wa(wa<0) = 0.01; 

            % ------------------ Call SiTHv2 ------------------------------
            [ETi, Tri, Esi, Eii, Esbi, SMi, RFi, GWi, snpx] = cal_SiTHv2(Rni,...
                Tai, Tasi, Precii, Pai, Gi, LAIii, Top, s_VODi, ...
                soilpar, pftpar, wa, zgw, snp, spinfg);
            % ------------------ Call SiTHv2 ------------------------------
            
            % writeout
            X_vals(:, j, 1) = ETi;
            X_vals(:, j, 2) = Tri;
            X_vals(:, j, 3) = Esi;
            X_vals(:, j, 4) = Eii;
            X_vals(:, j, 5) = Esbi;
            X_vals(:, j, 6) = SMi(:, 1);
            X_vals(:, j, 7) = SMi(:, 2);
            X_vals(:, j, 8) = SMi(:, 3);
            X_vals(:, j, 9) = RFi;
            X_vals(:, j, 10) = GWi;

            X_upti(1, j, 1) = SMi(end, 1);
            X_upti(1, j, 2) = SMi(end, 2);
            X_upti(1, j, 3) = SMi(end, 3);
            X_upti(1, j, 4) = GWi(end, 1);
            X_upti(1, j, 5) = snpx;

        end

        X_ET(i, :, :) = permute(X_vals(:, :, 1), [3, 2, 1]); % ET
        X_Tr(i, :, :) = permute(X_vals(:, :, 2), [3, 2, 1]); % Tr
        X_Es(i, :, :) = permute(X_vals(:, :, 3), [3, 2, 1]); % Es
        X_Ei(i, :, :) = permute(X_vals(:, :, 4), [3, 2, 1]); % Ei
        X_Esb(i, :, :) = permute(X_vals(:, :, 5), [3, 2, 1]); % Esb
        X_SM1(i, :, :) = permute(X_vals(:, :, 6), [3, 2, 1]); % SM1
        X_SM2(i, :, :) = permute(X_vals(:, :, 7), [3, 2, 1]); % SM2
        X_SM3(i, :, :) = permute(X_vals(:, :, 8), [3, 2, 1]); % SM3
        X_RF(i, :, :) = permute(X_vals(:, :, 9), [3, 2, 1]); % RF
        X_GW(i, :, :) = permute(X_vals(:, :, 10), [3, 2, 1]); % GW

        X_upt(i, :, :) = X_upti;

        ppm.increment();
    end

    % toc
    % Delete the progress handle when the parfor loop is done.
    delete(ppm);

    % save update variables to current folder
    save('upt_vals_01.mat', 'X_upt');

    % save the variables to output folder
    disp('------------------------- Writing results to mat files ... ')
    filename = [Output_folder 'SiTHv2.ET.S100.A' num2str(yr) '.mat'];
    X_ET = int16(100.*X_ET);
    save(filename, 'X_ET', '-v7.3');

    filename = [Output_folder 'SiTHv2.Tr.S100.A' num2str(yr) '.mat'];
    X_Tr = int16(100.*X_Tr);
    save(filename, 'X_Tr', '-v7.3');

    filename = [Output_folder 'SiTHv2.Es.S100.A' num2str(yr) '.mat'];
    X_Es = int16(100.*X_Es);
    save(filename, 'X_Es', '-v7.3');

    filename = [Output_folder 'SiTHv2.Ei.S100.A' num2str(yr) '.mat'];
    X_Ei = int16(100.*X_Ei);
    save(filename, 'X_Ei', '-v7.3');

    % filename = [Output_folder 'SiTHv2.En.S100.A' num2str(yr) '.mat'];
    % X_Esb = int16(100.*X_Esb);
    % save(filename, 'X_Esb', '-v7.3');

    % filename = [Output_folder 'SiTHv2.RF.S10.A' num2str(yr) '.mat']; % mm
    % X_RF = int16(10.*X_GW);
    % % imagesc(sum(X_RF,3))
    % save(filename, 'X_RF', '-v7.3');

    % filename = [Output_folder 'SiTHv2.GW.S100.A' num2str(yr) '.mat']; % mm to m
    % X_GW = int16(X_GW./10);
    % save(filename, 'X_GW', '-v7.3');

    filename = [Output_folder 'SiTHv2.SM.S100.A' num2str(yr) '.mat'];
    X_SM1 = int16(10000.*X_SM1); 
    X_SM2 = int16(10000.*X_SM2); 
    X_SM3 = int16(10000.*X_SM3); 
    save(filename, 'X_SM1', 'X_SM2', 'X_SM3', '-v7.3');

    disp('------------------------- End of this year ... ')
    disp(' ')

end
