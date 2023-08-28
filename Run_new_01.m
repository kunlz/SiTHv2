
clc
clear 
close all

res = 0.1;
r_m = 180 / res;
r_n = 360 / res;
% Output files folder :
Output_folder = 'Y:/SiTHv2_out_longterm/ERA5Lnew2/';

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
% subfolder_SM = 'SMcci/SM.CCIv071.scale00001.A'; 


subfolder_LAI = 'LAI/GEOV2_LAI/THEIA_GEOV2_R01_AVHRR_LAI_A'; 
subfolder_LC = 'LC/hildap_vGLOB1.0f/hildap_vGLOB.A';
% subfolder_VOD = 'Vegetation_Optical_Depth/VODCA/0.1deg/VODCA_Xband_A';

% subfolder_IWU = 'IWU/0.1deg/IWU_d_ens_scale001_A';

% georeference
latlim = [-90,90];
lonlim = [-180,180];
rasterSize = [720,1440];
RA = georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');

% Load the Soil Type dataset
Soilraster = load('inpara\Soilraster.mat');
Soilraster = Soilraster.Soilraster; 
% [Soilraster,~] = georesize(Soilraster,RA,2.5,"nearest"); 

% load the land mask
maskland = load('inpara\mask01new2.mat');
maskland = uint8(maskland.mask2);

% Load the optimal temperature for plant growth
Topt = load('inpara\Topt_v4.mat');
Topt = single(Topt.Topt_new);

% Parallel calculation initialization
% parpool('local', 40);

% Main loops
for yr = 1982 : 2020

    clear waa zgww snpp
    disp(' ')
    disp(['------------------------- Start calculation for ' num2str(yr)])
    % spin-up year : 100 years
    if yr == 1981

        spinfg = 1;
        disp(' ')
        disp('start year ... spin-up ... set spinfg = 1')

        % Initialization
        waa = 0.25 .* ones(r_m, r_n, 3); % initial value for swc
        zgww = 5050 .* ones(r_m, r_n); % initial value for groundwater table
        snpp = zeros(r_m, r_n); % initial value for snowpack depth

        % load the updated variables
%         uptval = load('upt_vals_01_backup.mat');
%         uptval = uptval.X_upt;
%         
%         waa = uptval(:,:,1:3); 
%         zgww = uptval(:,:,4); 
%         snpp = uptval(:,:,5); 

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
    % load forcing data %
    % ----------------- %
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
%     disp('Construct Object for Satellite-based VOD  ...')
%     EMO_VOD = matfile([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    % EMO_VOD = load([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    % EMO_VOD = EMO_VOD.VODCAy;

    % Satellite-based Landcover
    disp(['Load Satellite-based Landcover  ... For the year :: ' num2str(yr)])
    LC_year = load([Forcing_folder subfolder_LC num2str(yr) '.mat']);
    LC_year = LC_year.LULC; 

    % load IWU
    % EMO_IWU = matfile([Forcing_folder subfolder_IWU num2str(yr) '.mat']);
    % IWU = single(IWU.IWUx);

    % [~, ~, days] = size(EMO_Rn);
    days = yeardays(yr);

    % ------------------ %
    % Parallel Computing %
    % ------------------ %

    disp('Preallocate memory to each variables ... ')
    
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
        
        % matfile read each row, latitude read
        Rnix = permute(EMO_Rn(i, :, :), [3, 2, 1]);

        Taix = permute(EMO_Ta(i, :, :), [3, 2, 1]);

        Precix = permute(EMO_Preci(i, :, :), [3, 2, 1]);

        Paix = permute(EMO_Pa(i, :, :), [3, 2, 1]);

        LAIix = permute(EMO_LAI(i, :, :), [3, 2, 1]);

        % VODix = permute(EMO_VOD(i, :, :), [3, 2, 1]);

        % read IWU
        % IWUix = permute(EMO_IWU.IWUx(i, :, :), [3, 2, 1]);

        % X_mat file for different output variables
        X_vals = zeros(days, r_n, 10);
        X_upti = zeros(1, r_n, 5);

        for j = 1 : r_n

            if maskland(i, j) == 0
                continue
            end

            % update variables
            wa = reshape(waa(i, j, :), [1, 3]);
            zgw = zgww(i, j);
            snp = snpp(i, j);

            % Meteo forcing for each pixel, rescale
            Rni = 0.01.*double(Rnix(:, j));
            Tai = 0.01.*double(Taix(:, j));
            Precii = 0.01.*double(Precix(:, j));
            % IWUi = 0.01.*double(IWUix(:,j));
            IWUi = zeros(days,1);
            Pai = 0.01.*double(Paix(:, j));

            % calculate Tas
            Tasi = Tai;
            Tasi(Tasi < 0) = 0;
            Tasi = cumsum(Tasi);

            % Satellite LAI for each pixel
            LAIi = 0.01.*double(LAIix(:, j));
            xo = day(EMO_LAItime,"dayofyear");
            xi = 1:1:days;

            LAIii = interp1(xo', LAIi, xi', 'pchip', 'extrap');
            % plot(xo',LAIi,'o',xi',LAIii);
            LAIii(LAIii < 0) = 0.2;

            % Cal G_soil, % Choudhury et al., 1987
            Gi = 0.4 .* Rni .* exp(-0.5 .* LAIii);

            % Get Stress from VOD
%             VODi = 0.001.*double(VODix(1:days, j));
%             VODi = smooth(VODi, 7, 'moving');
%             VODi(VODi < 0) = 0;
%             s_VODi = (VODi ./ max(VODi)).^0.5;
            s_VODi = ones(days,1);

            % Get Topt
            Top = Topt(i, j); % Top = 25;
            if isnan(Top)
                Top = 25;
            end

            % Parameters set for plant and soil
            % 1- PlantType
            PFTi = LC_year(i, j);
            pftpar = get_pftpar_new(PFTi);

            % 2- SoilType
            SC = Soilraster(i, j);
            if SC == 0
                SC = 6;
            end
            soilpar = get_soilpar_raster(SC);

            % optpara
            optpara = get_optpara_new(PFTi);
            optpara(1) = 1.26;
            
            % check wa
            wa(wa<0) = 0.01; 

            % ------------------ Call SiTH --------------------------------
            [ETi, Tri, Esi, Eii, Esbi, SMi, RFi, GWi, snpx] = cal_SiTH(Rni,...
                Tai, Tasi, Precii, IWUi, Pai, Gi, LAIii, Top, s_VODi, ...
                soilpar, pftpar, wa, zgw, snp, optpara, spinfg);
            % ------------------ Call SiTH --------------------------------
            
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
    filename = [Output_folder 'SiTH.ET.S100.A' num2str(yr) '.mat'];
    X_ET = int16(100.*X_ET);
    save(filename, 'X_ET', '-v7.3');

    filename = [Output_folder 'SiTH.Tr.S100.A' num2str(yr) '.mat'];
    X_Tr = int16(100.*X_Tr);
    save(filename, 'X_Tr', '-v7.3');

    filename = [Output_folder 'SiTH.Es.S100.A' num2str(yr) '.mat'];
    X_Es = int16(100.*X_Es);
    save(filename, 'X_Es', '-v7.3');

    filename = [Output_folder 'SiTH.Ei.S100.A' num2str(yr) '.mat'];
    X_Ei = int16(100.*X_Ei);
    save(filename, 'X_Ei', '-v7.3');

    filename = [Output_folder 'SiTH.En.S100.A' num2str(yr) '.mat'];
    X_Esb = int16(100.*X_Esb);
    save(filename, 'X_Esb', '-v7.3');

%     filename = [Output_folder 'SiTH.RF.S10.A' num2str(yr) '.mat']; % mm
%     X_RF = int16(10.*X_GW);
%     % imagesc(sum(X_RF,3))
%     save(filename, 'X_RF', '-v7.3');

%     filename = [Output_folder 'SiTH.GW.S100.A' num2str(yr) '.mat']; % mm to m
%     X_GW = int16(X_GW./10);
%     save(filename, 'X_GW', '-v7.3');

    filename = [Output_folder 'SiTH.SM.S100.A' num2str(yr) '.mat'];
    X_SM1 = int16(10000.*X_SM1); 
    X_SM2 = int16(10000.*X_SM2); 
    X_SM3 = int16(10000.*X_SM3); 
    save(filename, 'X_SM1', 'X_SM2', 'X_SM3', '-v7.3');

    % writeout to nc files
    %     disp('------------------------- Writing results to netCDF files ... ')
    %     filename = [Output_folder 'SiTH.ET.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'Et', 'Actual Evapotranspiration', X_ET, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.Tr.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'Tr', 'Plant Transpiration', X_Tr, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.Es.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'Es', 'Soil Evaporation', X_Es, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.Ei.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'Ei', 'Interception Evaporation', X_Ei, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.En.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'En', 'Snow\Ice Sublimation', X_Esb, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.RF.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'RF', 'Surface Runoff', X_RF, days, 'mm');
    %
    %     filename = [Output_folder 'SiTH.GW.A' num2str(yr) '.nc'];
    %     write2nc(filename, 'GW', 'Groundwater table depth', X_GW, days, 'm');
    %     tic
    %     filename = [Output_folder 'SiTH.SM.A' num2str(yr) '.nc'];
    %     write2nc3SM(filename, 'sm1', 'Soil Moisture at layer #1', ...
    %         'sm2', 'Soil Moisture at layer #2', 'sm3', 'Soil Moisture at layer #3', ...
    %         X_SM1, X_SM2, X_SM3, days, '0.01');
    %     toc

    disp('------------------------- End of this year ... ')
    disp(' ')

end
