

load testdata_arou_v1.mat

Rain    =   data(:,1);
ETo     =   data(:,2);
Ta      =   data(:,3);
Pa      =   data(:,4);
Rn      =   data(:,5);
G       =   data(:,6);
SWC1    =   data(:,7);
SWC2    =   data(:,8);
LAI     =   data(:,9);


% NDVI = data(:,end);
% fc = NDVI-0.05;               % Fraction total vegetation cover
% LAI = -log(1-fc)./0.5; 


soilpar = [0.5353, -0.0316, 0.416, 4.5, 0.265, 0.208, 0.087];
pftpar = [0.01, 13, -1.92, 2.6];
wa = [0.5, 0.5, 0.5];
zgw = 20;
SnpN = 2;

% for i = 1 : 365
%     i
%     
%     Rni = Rn(i,1);
%     Tai = Ta(i,1);
%     Raini = Rain(i,1);
%     Pai = Pa(i,1);
%     Gi = G(i,1);
%     LAIi = LAI(i,1);
%     
% 
%     [Et, Tr, Es, Ei, Esnow, wa, srf, zgw, SnpN] = STHM(Rni, Tai, 25, Raini, Pai, Gi, LAIi, soilpar, pftpar, wa, zgw, SnpN);
%     
% end

for i = 1 : 92
    i
    
    Rni = Rn(i,1);
    Tai = Ta(i,1);
    Raini = Rain(i,1);
    Pai = Pa(i,1);
    Gi = G(i,1);
    LAIi = LAI(i,1);
    

%     [Et(i,1), Tr(i,1), Es(i,1), Ei(i,1), Esnow(i,1), wa, srf, zgw, SnpN] = STHM(Rni, Tai, 25, Raini, Pai, Gi, LAIi, soilpar, pftpar, wa, zgw, SnpN);
    
        % Set the soil depth for three soil layers
    zm = [50, 500, 2000]; % mm

    % Potential Evaporation allocated to canopy and soil surface
    [pEc(i,1), pEs(i,1)] = potentialET(Rni, Gi, LAIi, Tai, Pai);

    % Interception evaporation
    [Ei, wet] = interception(LAIi, Raini, pEc(i,1), pftpar);

    % Snow sublimation
    [SnpN, Esnow, snowmelt] = snowpack_balance(Raini, Tai, SnpN);

    % Calculate the net Precipitation into soil surface
    Pnet = max(0, Raini + snowmelt - Esnow - Ei);

    % Calculate the Runoff
    [srf, IWS] = runoff(Pnet, zgw, zm, wa, soilpar);

    % Calculate Plant transpiration, Soil evaporation, and Soil moisture status
    [wa, zgw, Tr, Es, uex] = sw_balance(IWS, pEc(i,1), pEs(i,1), Tai, 25, wa, soilpar, pftpar, wet, zm, zgw);

    % Calculate the total Evaporanspiration
    Et(i,1) = Tr + Es + Ei + Esnow;
    srf = srf + uex;
    
    wo(i,:) = wa;
    zzg(i,1) = zgw/1000;
    
end

plot([ETo,Et]);

scatter(ETo,Et);
axis square
xlim([0,5]);
ylim([0,5]); 

close all
plot([SWC1]);
yyaxis right
plot(wo(:,3))

corr(SWC1, wo(:,2))

corr(ETo,Et)










