close all
clearvars
clc

fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\9. OSM OFDM\';
        
end
figname = 'SIM_03_SNR vs Offset';
% fname = [ctDirRes figname];
fname = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\9. OSM OFDM\SIM_03_SNR vs Offset';
% Y = 1:2:20;
% X = 1:numel(Y);
uiopen([fname '.fig'],1);

fig = gcf;
title('SNR vs Offset, Target BER = 10^{-3}');
set(fig,'DefaultLineLineWidth',2);
set(fig,'DefaultAxesFontName','Helvetica');
set(fig,'DefaultAxesFontSize',16);
% strLgd = [{'SM: \mu_s=0.50'}...
%     {'SMP: \mu_s=0.00'},...
%     {'SM: \mu_s=1.00'},...
%     {'SMP: \mu_s=1.00'},...
%     {'SM: \mu_s=1.41'},...
%     {'SMP: \mu_s=1.41'},...
%     {'SM: \mu_s=2.00'},...
%     {'SMP: \mu_s=2.00'}];
% legend(strLgd,'Location','SouthWest');
% plot(X,Y);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 4]);

saveas(fig,[fname '.fig'],'fig');
saveas(fig,[fname '.png'],'png');
% saveas(fig,[fname '.eps'],'eps');