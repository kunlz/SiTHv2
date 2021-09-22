clear
clc
close all

% Output files folder :
Output_folder = '../SiTH_out/';

% Forcing data folder :
Focing_folder = '../SiTHM_ForcingData/';

% Sub-folders of Focing data
subfolder_Rn = 'Net_Radiation/ERA5_land_surface_net_solar_radiation_';
subfolder_Ta = 'Air_Temperature/ERA5_2m_temperature_';
subfolder_Pa = 'Air_Pressure/ERA5_land_surface_pressure_';

subfolder_LAI = 'LAI/GLASS01B02.V40.A';
subfolder_LC = 'Landcover/MCD12C1.A';
subfolder_VOD = 'Vegetation_Optical_Depth/VODCA/VODCA_Xband_A';
subfolder_Preci = 'Precipitation/MSWEP/MSWEP_V280_A';

% Load the Soil Type dataset
Soilraster = load([Focing_folder 'Soilraster.mat']);
Soilraster = Soilraster.Soilraster;

% load the land mask
maskland = load([Focing_folder 'maskland.mat']);
maskland = maskland.maskland;

% Load the optimal temperature for plant growth
Topt = load([Focing_folder 'Topt_v2.mat']);
Topt = Topt.Topt;

% Parallel calculation initialization
% parpool('local', 10);

% Main loops
for yr = 2001:2018
    
    clear waa zgww snpp
    disp(' ')
    disp(['------------------------- Start calculation for ' num2str(yr)] )
    % spin-up year : 20 years
    if yr == 2001

        spinfg = 0;
        disp(' ')
        disp('normal year ... set spinfg = 0')
        
        % Initialization
%         waa = 0.3 .* ones(1800, 3600, 3); % initial value for swc
%         zgww = 5050 .* ones(1800, 3600); % initial value for groundwater table
%         snpp = zeros(1800, 3600); % initial value for snowpack depth

        % load the updated variables
        uptval = load('upt_vals_500yr.mat');
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
        
        waa = uptval(:,:,1:3);
        zgww = uptval(:,:,4);
        snpp = uptval(:,:,5);
    end

    % ----------------- %
    % load forcing data %
    % ----------------- %
    % Net Radiation, W/m2
    disp('Construct Object for Net Radiation ...')
    EMO_Rn = matfile([Focing_folder subfolder_Rn num2str(yr) '.mat']);

    % Air Temperature, C, 2m
    disp('Construct Object for Air Temperature ...')
    EMO_Ta = matfile([Focing_folder subfolder_Ta num2str(yr) '.mat']);

    % Precipitation, mm
    disp('Construct Object for Precipitation ...')
    EMO_Preci = matfile([Focing_folder subfolder_Preci num2str(yr) '.mat']);

    % Air Pressure, kPa
    disp('Construct Object for Air Pressure ...')
    EMO_Pa = matfile([Focing_folder subfolder_Pa num2str(yr) '.mat']);

    % Satellite-based LAI
    disp('Construct Object for Satellite-based LAI  ...')
    EMO_LAI = matfile([Focing_folder subfolder_LAI num2str(yr) '.mat']);

    % Satellite-based VOD
    disp('Construct Object for Satellite-based VOD  ...')
    EMO_VOD = matfile([Focing_folder subfolder_VOD num2str(yr) '.mat']);

    % Satellite-based Landcover
    disp('Load Satellite-based Landcover  ...')
    LC_year = load([Focing_folder subfolder_LC num2str(yr) '001.006.mat']);
    LC_year = LC_year.LC_IGBP;

    [~, ~, days] = size(EMO_Rn, 'ERA5L_Rn');

    % ------------------ %
    % Parallel Computing %
    % ------------------ %
    
    disp('Preallocate memory to each variables ... ')
%     X_ET = int16(zeros(1800, 3600, days));
%     X_Tr = int16(zeros(1800, 3600, days));
%     X_Es = int16(zeros(1800, 3600, days));
%     X_Ei = int16(zeros(1800, 3600, days));
%     X_Esb = int16(zeros(1800, 3600, days));
    X_SM1 = int8(zeros(1800, 3600, days));
    X_SM2 = int8(zeros(1800, 3600, days));
    X_SM3 = int8(zeros(1800, 3600, days));
    X_RF = int16(zeros(1800, 3600, days));
    X_GW = int16(zeros(1800, 3600, days));
    
    X_upt = zeros(1800, 3600,5);

    ppm = ParforProgressbar(1800, 'showWorkerProgress', true);
    % tic
    parfor i = 1:1800

        % matfile read each row, latitude read
        Rnix = permute(EMO_Rn.ERA5L_Rn(i, :, 1:days), [3, 2, 1]);
        Rnix = 0.01 .* double(Rnix);

        Taix = permute(EMO_Ta.ERA5L_Ta(i, :, 1:days), [3, 2, 1]);
        Taix = 0.01 .* double(Taix);

        Precix = permute(EMO_Preci.MSWEP(i, :, 1:days), [3, 2, 1]);
        Precix = 0.01 .* double(Precix);

        Paix = permute(EMO_Pa.ERA5L_Pa(i, :, 1:days), [3, 2, 1]);
        Paix = 0.01 .* double(Paix);

        LAIix = permute(EMO_LAI.GLASSLAI(i, :, :), [3, 2, 1]);
        LAIix = 0.01 .* double(LAIix);

        VODix = permute(EMO_VOD.VODCAy(i, :, 1:days), [3, 2, 1]);
        VODix = 0.001 .* double(VODix);
        
        % X_mat file for different output variables
        X_vals = zeros(days, 3600, 5);
        X_upti = zeros(1, 3600, 5);

        for j = 1:3600

            if LC_year(i, j) == 0
                continue
            end

            % update variables
            wa = reshape(waa(i, j, :), [1, 3]);
            zgw = zgww(i, j);
            snp = snpp(i, j);

            % Meteo forcing for each pixel
            Rni = Rnix(:, j);
            Tai = Taix(:, j);
            Precii = Precix(:, j);
            Pai = Paix(:, j);

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
            VODi = VODix(:, j);
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

            % Call the SiTHM
            [ETi, Tri, Esi, Eii, Esbi, SMi, RFi, GWi, snpx] = Cal_SiTHM(Rni, Tai, ...
                Precii, Pai, Gi, LAIii, Top, s_VODi, ...
                soilpar, pftpar, wa, zgw, snp, optpara, spinfg);

            % writeout
%             X_vals(:, j, 1) = ETi;
%             X_vals(:, j, 2) = Tri;
%             X_vals(:, j, 3) = Esi;
%             X_vals(:, j, 4) = Eii;
%             X_vals(:, j, 5) = Esbi;
            X_vals(:, j, 1) = SMi(:, 1);
            X_vals(:, j, 2) = SMi(:, 2);
            X_vals(:, j, 3) = SMi(:, 3);
            X_vals(:, j, 4) = RFi;
            X_vals(:, j, 5) = GWi;
            

            X_upti(1,j,1) = SMi(end,1);
            X_upti(1,j,2) = SMi(end,2);
            X_upti(1,j,3) = SMi(end,3);
            X_upti(1,j,4) = GWi(end,1);
            X_upti(1,j,5) = snpx;

        end

%         X_ET(i, :, :) = permute(int16(100 .* X_vals(:, :, 1)), [3, 2, 1]); % ET  --> /100 (12.34 mm) 
%         X_Tr(i, :, :) = permute(int16(100 .* X_vals(:, :, 2)), [3, 2, 1]); % Tr  --> /100 (12.34 mm) 
%         X_Es(i, :, :) = permute(int16(100 .* X_vals(:, :, 3)), [3, 2, 1]); % Es  --> /100 (12.34 mm) 
%         X_Ei(i, :, :) = permute(int16(100 .* X_vals(:, :, 4)), [3, 2, 1]); % Ei  --> /100 (12.34 mm) 
%         X_Esb(i, :, :) = permute(int16(100 .* X_vals(:, :, 5)), [3, 2, 1]); % Esb  --> /100 (12.34 mm) 
        X_SM1(i, :, :) = permute(int8(100 .* X_vals(:, :, 1)), [3, 2, 1]); % SM1  --> /100  (0.12)
        X_SM2(i, :, :) = permute(int8(100 .* X_vals(:, :, 2)), [3, 2, 1]); % SM2  --> /100  (0.12)
        X_SM3(i, :, :) = permute(int8(100 .* X_vals(:, :, 3)), [3, 2, 1]); % SM3  --> /100  (0.12)
        X_RF(i, :, :) = permute(int16(X_vals(:, :, 4)), [3, 2, 1]); % RF,  --> /1 (12 mm) 
        X_GW(i, :, :) = permute(int16(0.1 .* X_vals(:, :, 5)), [3, 2, 1]); % GW ,  --> /100 (12.34 m) 
        
        X_upt(i, :, :) = X_upti;
        
        % Writeout to file
        % rows = num2str(i,'%04d');
        % outname = [Output_folder '\tem.A' num2str(yr) '.R' rows '.mat'];
        % parsave(outname, X_vals);
        ppm.increment();
    end
    % toc
    % Delete the progress handle when the parfor loop is done.
    delete(ppm);
    
    % save update variables to current folder
    save('upt_vals.mat', 'X_upt');
    
    % save the variables to output folder
    disp('------------------------- Writing results to mat files ... ')
%     filename = [Output_folder 'SiTH.ET.A' num2str(yr) '.mat'];
%     save(filename, 'X_ET', '-v7.3');
%     
%     filename = [Output_folder 'SiTH.Tr.A' num2str(yr) '.mat'];
%     save(filename, 'X_Tr', '-v7.3');
%     
%     filename = [Output_folder 'SiTH.Es.A' num2str(yr) '.mat'];
%     save(filename, 'X_Es', '-v7.3');
%     
%     filename = [Output_folder 'SiTH.Ei.A' num2str(yr) '.mat'];
%     save(filename, 'X_Ei', '-v7.3');
%     
%     filename = [Output_folder 'SiTH.En.A' num2str(yr) '.mat'];
%     save(filename, 'X_Esb', '-v7.3');
    
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

%             figure(1)
%             plot([ETi, Tri, Esi, Eii, Esbi])
%             figure(2)
%             plot(SMi)
%             figure(3)
%             plot(RFi)
%             figure(4)
%             plot(GWi)

% ETs = sum(0.01.*double(X_ET),3);
%
% figure(6)
% imagesc(ETs)

% figure(2)
% plot(VPDi)
%
% figure(3)
% plot(ETi)
% sum(ETi)
%
%
% figure(4)
% plot(ETi)
% sum(ETi)
