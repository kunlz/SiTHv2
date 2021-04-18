
clear
clc

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

soilpar = get_soilpar(5);
pftpar = get_pftpar(10); % grassland

wa = [0.3, 0.3, 0.3]; % initial value for swc
zgw = 200; % initial value for groundwater table
snp = 0; % initial value for snowpack depth
%%
% set spin-up times
m = 20;

for nn = 1 : m
    for i = 1 : 92
        
        Rni = Rn(i,1);
        Tai = Ta(i,1);
        Pei = Rain(i,1);
        Pai = Pa(i,1);
        Gi = G(i,1);
        LAIi = LAI(i,1);
        
        [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp] = SiTHM(Rni, Tai, 25,...
            Pei, Pai, Gi, LAIi, soilpar, pftpar, wa, zgw, snp);
        
        
        ETs(i,1) = Et;
        wo(i,:) = wa;
        zzg(i,1) = zgw/1000;
        runof(i,1) = srf;
        
    end
end
%% Plot section

close all
set(gcf,'unit','centimeters','position',[10,1,30,18]);
ha = tight_subplot(2, 2, [0.1 0.1], [.1 .1], [.08 .08]);

axes(ha(1));
plot(ETo, 'DisplayName','Observed ET','LineWidth',1.5);
hold on
plot(ETs, 'DisplayName','Simulated ET','LineWidth',1.5);
ylabel('ET (mm day^{-1})')
legend({'Observed ET','Simulated ET'});
set(gca,'FontSize',12);
title('Evapotranspiration')

axes(ha(2));
plot(SWC1,'DisplayName','Observed SWC','LineWidth',1.5,...
    'Color', [0.47,0.67,0.19]);
set(gca,'Ycolor',[0.47,0.67,0.19])
ylabel('Observed SWC (%)')

yyaxis right
plot(wo(:,3) .* 100,'DisplayName','Simulated SWC','LineWidth',1.5,...
    'Color', [0.49,0.18,0.56]);
set(gca,'Ycolor',[0.49,0.18,0.56])
ylabel('Simulated SWC (%)','Rotation',-90)

set(gca,'FontSize',12);
title('Soil Water Content')

legend({'Observed SWC','Simulated SWC'});

axes(ha(3));
plot(zzg,'DisplayName','Groundwater Table Depth','LineWidth',1,...
    'Color', [0.15,0.15,0.15]);
ylim([0,3]);
ylabel('Groundwater Table Depth (m)')
set(gca,'Ycolor',[0.15,0.15,0.15],'FontSize',12)
legend({'Groundwater Table Depth'});
title('Groundwater')

axes(ha(4));
plot(runof,'DisplayName','Surface runoff','LineWidth',1.5,...
    'Color', [0.64,0.08,0.18])
ylabel('Surface runoff (mm)')
set(gca,'Ycolor',[0.64,0.08,0.18])
yyaxis right
bar(Rain,'DisplayName','Rainfall','LineWidth',1.5,...
    'FaceColor', [0.00,0.45,0.74], 'EdgeColor','none', 'FaceAlpha',0.6);
set(gca,'YDir','reverse','FontSize',12);
set(gca,'Ycolor',[0.00,0.45,0.74])
ylabel('Rainfall (mm)','Rotation',270)
legend({'Surface runoff','Rainfall'});
title('Runoff & Rainfall')

%%
% export figure
exportgraphics(gcf,'testFigure_Arou.png','Resolution',600);










