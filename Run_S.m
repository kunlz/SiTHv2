clear
clc
close all

res = 0.25;

% Output files folder :
Output_folder = '../SiTH_out/';

if exist(Output_folder, 'dir') == 0.
    mkdir(Output_folder);
end

% Forcing data folder :
Forcing_folder = '../SiTH_ForcingData/';

% Sub-folders of Focing data
subfolder_Rn = 'Net_Radiation/CERES/0.25deg/CERES_RN_A';
subfolder_Ta = 'Air_Temperature/0.25deg/TA_ERA5L_A';
subfolder_Pa = 'Air_Pressure/0.25deg/PA_ERA5L_A';

subfolder_LAI = 'LAI/GLASS/0.25deg/LAI_GLASS_A';
subfolder_LC = 'Landcover/0.25deg/LC_IGBP_A';
subfolder_VOD = 'Vegetation_Optical_Depth/VODCA/0.25deg/VODCA_A';
subfolder_Preci = 'Precipitation/ERA5L/P_ERA5L_tot_A';

% Load the Soil Type data
Soilraster = load([Forcing_folder 'Soilraster.mat']);
Soilraster = Soilraster.x2;

% Load the alfset data
alfset = load([Forcing_folder 'alfsets.mat']);
alfset = alfset.alfset;

% load the land mask
maskland = load([Forcing_folder 'maskland25.mat']);
maskland = single(maskland.mask25);

% Load the optimal temperature for plant growth
Topt = load([Forcing_folder 'Topt25_v2.mat']);
Topt = single(Topt.Topt25);

% Parallel calculation initialization
% parpool('local', 6);

% Main loops
for yr = 2001 : 2018

    clear waa zgww snpp
    disp(' ')
    disp(['------------------------- Start calculation for ' num2str(yr)])
    % spin-up year : 50 years
    if yr == 2000

        spinfg = 1;
        disp(' ')
        disp('start year ... set spinfg = 1')

        % Initialization
%         waa = 0.25 .* ones(720, 1440, 3); % initial value for swc
%         zgww = 5050 .* ones(720, 1440); % initial value for groundwater table
%         snpp = zeros(720, 1440); % initial value for snowpack depth

        % load the updated variables
        uptval = load('upt_vals.mat');
        uptval = uptval.X_upt;
        
        waa = uptval(:,:,1:3);
        zgww = uptval(:,:,4);
        snpp = uptval(:,:,5);

    else
        spinfg = 0;
        disp(' ')
        disp('normal year ... set spinfg = 0')
        % load the updated variables
        uptval = load('upt_vals.mat');
        uptval = uptval.X_upt;

        waa = uptval(:, :, 1:3);
        zgww = uptval(:, :, 4);
        snpp = uptval(:, :, 5);
    end

    % ----------------- %
    % load forcing data %
    % ----------------- %
    % Net Radiation, W/m2
    disp('Construct Object for Net Radiation ...')
    % EMO_Rn = matfile([Forcing_folder subfolder_Rn num2str(yr) '.mat']);
    EMO_Rn = load([Forcing_folder subfolder_Rn num2str(yr) '.mat']);
    EMO_Rn = EMO_Rn.CERES_RN; 

    % Air Temperature, C, 2m
    disp('Construct Object for Air Temperature ...')
    % EMO_Ta = matfile([Forcing_folder subfolder_Ta num2str(yr) '.mat']);
    EMO_Ta = load([Forcing_folder subfolder_Ta num2str(yr) '.mat']);
    EMO_Ta = EMO_Ta.TA_ERA5L; 

    % Precipitation, mm
    disp('Construct Object for Precipitation ...')
    % EMO_Preci = matfile([Forcing_folder subfolder_Preci num2str(yr) '.mat']);
    EMO_Preci = load([Forcing_folder subfolder_Preci num2str(yr) '.mat']);
    EMO_Preci = EMO_Preci.P_ERA5L;

    % Air Pressure, kPa
    disp('Construct Object for Air Pressure ...')
    % EMO_Pa = matfile([Forcing_folder subfolder_Pa num2str(yr) '.mat']);
    EMO_Pa = load([Forcing_folder subfolder_Pa num2str(yr) '.mat']);
    EMO_Pa = EMO_Pa.PA_ERA5L;

    % Satellite-based LAI
    disp('Construct Object for Satellite-based LAI  ...')
    % EMO_LAI = matfile([Forcing_folder subfolder_LAI num2str(yr) '.mat']);
    EMO_LAI = load([Forcing_folder subfolder_LAI num2str(yr) '.mat']);
    EMO_LAI = EMO_LAI.LAI_GLASS;

    % Satellite-based VOD
    disp('Construct Object for Satellite-based VOD  ...')
    % EMO_VOD = matfile([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    EMO_VOD = load([Forcing_folder subfolder_VOD num2str(yr) '.mat']);
    EMO_VOD = EMO_VOD.VODCA;

    % Satellite-based Landcover
    disp('Load Satellite-based Landcover  ...')
    LC_year = load([Forcing_folder subfolder_LC num2str(yr) '.mat']);
    LC_year = LC_year.LC_IGBP;

    % load IWU
    IWU = load(['../SiTH_ForcingData/IWU/IWU_d_ens_A' num2str(yr) '.mat']);
    IWU = single(IWU.IWUd);

    [~, ~, days] = size(EMO_Rn); 

    % ------------------ %
    % Parallel Computing %
    % ------------------ %

    disp('Preallocate memory to each variables ... ')
    r_m   = 180 / res;
    r_n   = 360 / res;
    X_ET  = single(zeros(r_m, r_n, days));
    X_Tr  = single(zeros(r_m, r_n, days));
    X_Es  = single(zeros(r_m, r_n, days));
    X_Ei  = single(zeros(r_m, r_n, days));
    X_Esb = single(zeros(r_m, r_n, days));
    X_SM1 = single(zeros(r_m, r_n, days));
    X_SM2 = single(zeros(r_m, r_n, days));
    X_SM3 = single(zeros(r_m, r_n, days));
    X_RF  = single(zeros(r_m, r_n, days));
    X_GW  = single(zeros(r_m, r_n, days));
    
    disp('Start calculation ... ')

    X_upt = zeros(r_m, r_n, 5);

    ppm = ParforProgressbar(r_m, 'showWorkerProgress', true);
    parfor i = 1:r_m

        % matfile read each row, latitude read
        Rnix = permute(EMO_Rn(i, :, :), [3, 2, 1]);

        Taix = permute(EMO_Ta(i, :, :), [3, 2, 1]);

        Precix = permute(EMO_Preci(i, :, :), [3, 2, 1]);

        Paix = permute(EMO_Pa(i, :, :), [3, 2, 1]);

        LAIix = permute(EMO_LAI(i, :, :), [3, 2, 1]);

        VODix = permute(EMO_VOD(i, :, :), [3, 2, 1]);

        % read IWU
        IWUix = permute(IWU(i, :, :), [3, 2, 1]);

        % X_mat file for different output variables
        X_vals = zeros(days, r_n, 10);
        X_upti = zeros(1, r_n, 5);

        for j = 1:r_n

            if LC_year(i, j) == 0 || maskland(i, j) == 0
                continue
            end

            % update variables
            wa = reshape(waa(i, j, :), [1, 3]);
            zgw = zgww(i, j);
            snp = snpp(i, j);

            % Meteo forcing for each pixel
            Rni = Rnix(1:days, j);
            Tai = Taix(1:days, j);
            Precii = Precix(1:days, j);
            IWUi = IWUix(1:days,j);
            % IWUi = zeros(length(IWUi),1);
            Pai = Paix(1:days, j);

            % calculate Tas
            Tasi = Tai;
            Tasi(Tasi < 0) = 0;
            Tasi = cumsum(Tasi);

            % Satellite LAI for each pixel
            LAIi = LAIix(:, j);
            xo = 1:8:366;
            xi = 1:1:days;
            LAIii = interp1(xo', LAIi, xi', 'pchip', 'extrap');
            % plot(xo',LAIi,'o',xi',LAIii);
            LAIii(LAIii < 0) = 0.2;

            % Cal G_soil, % Choudhury et al., 1987
            Gi = 0.4 .* Rni .* exp(-0.5 .* LAIii);

            % Get Stress from VOD
            VODi = VODix(1:days, j);
            VODi = smooth(VODi, 7, 'moving');
            VODi(VODi < 0) = 0;
            s_VODi = (VODi ./ max(VODi)).^0.5;

            % Get Topt
            Top = Topt(i, j); % Top = 25;

            % Parameters set for plant and soil
            % 1- PlantType
            PFTi = LC_year(i, j);
            pftpar = get_pftpar_raster(PFTi);

            % 2- SoilType
            SC = Soilraster(i, j);
            soilpar = get_soilpar_raster(SC);

            % optpara
            optpara = get_optpara(PFTi);
            alf = alfset(i,j);
            if isnan(alf)
                alf = 1.26;
            end
            optpara(1) = alf;
            
            % check wa
            wwp = soilpar(7);
            wa(wa<wwp) = wwp +0.001; 

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
    save('upt_vals.mat', 'X_upt');

    % save the variables to output folder
    disp('------------------------- Writing results to mat files ... ')
    filename = [Output_folder 'SiTH.ET.A' num2str(yr) '.mat'];
    save(filename, 'X_ET', '-v7.3');

    filename = [Output_folder 'SiTH.Tr.A' num2str(yr) '.mat'];
    save(filename, 'X_Tr', '-v7.3');

    filename = [Output_folder 'SiTH.Es.A' num2str(yr) '.mat'];
    save(filename, 'X_Es', '-v7.3');

    filename = [Output_folder 'SiTH.Ei.A' num2str(yr) '.mat'];
    save(filename, 'X_Ei', '-v7.3');

    filename = [Output_folder 'SiTH.En.A' num2str(yr) '.mat'];
    save(filename, 'X_Esb', '-v7.3');

    filename = [Output_folder 'SiTH.RF.A' num2str(yr) '.mat'];
    save(filename, 'X_RF', '-v7.3');

    filename = [Output_folder 'SiTH.GW.A' num2str(yr) '.mat'];
    save(filename, 'X_GW', '-v7.3');

    filename = [Output_folder 'SiTH.SM.A' num2str(yr) '.mat'];
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
