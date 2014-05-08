% scrOFDMWDMPL
% script to plot data generated by scrOFDMWDMPL
function scrOFDMWDMPL(dataFile)
if ~exist('dataFile','var')
    error('Data file not specified');
end
close all;
clearvars -except dataFile;
clc;
load(dataFile);

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
% daxesfontsize = get(0,'DefaultAxesFontSize');
% set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');
dfigppm = get(0,'DefaultFigurePaperPositionMode');
set(0,'DefaultFigurePaperPositionMode','Manual');
dfigpu = get(0,'DefaultFigurePaperUnits');
set(0,'DefaultFigurePaperUnits','Inches');
dfigpp = get(0,'DefaultFigurePaperPosition');
set(0,'DefaultFigurePaperPosition',[0 0 11 8.5]);
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',6);
FIGTITLE = 'Off';

h = waitbar(0,'Plotting Results: 0.00% done','Name',WBTITLE,'Visible','Off');
set(h,'Position',[WBX WBY WBW WBH],'Visible','On');

%% PLOT
% PLOT Configs
PLOTDMIN = 4;
PLACOLC = 'm'; PLACOLS = '-'; PLACOMK = 'o';
PLDCOLC = 'c'; PLDCOLS = '-'; PLDCOMK = 'd';
PLTXLCS = {'r';'g';'b'}; PLTXLSS = {'--';'-.';':'}; PLTXMKS = {'>';'s';'*'};
% Figure CCT config
FIGCCT = figure('Name',sprintf('PSD vs CCT'),'NumberTitle',FIGTITLE);
FIGCCTPLNR = power(2,floor(log2(LENCCTPL)/2));
FIGCCTPLNC = power(2,ceil(log2(LENCCTPL)/2));
FIGLDAMIN = 400; FIGLDAMAX = 800;
iTPL = 1;
% Figure BER vs SNR config
FIGBER = zeros(1,LENCCT);
FIGBERPLNR = power(2,ceil(log2(LENOFDMTYPES)/2));
FIGBERPLNC = power(2,floor(log2(LENOFDMTYPES)/2));
FIGBERXMIN = min(RNGSNRDB); FIGBERXMAX = max(RNGSNRDB);
FIGBERYMIN = 0.9*BERTH; FIGBERYMAX = 1;
FIGBERLGD = {};

TOTALLOOPS = LENCCT*LENOFDMTYPES+3;
LOOPCOUNT = 0;
for iT = 1:LENCCT
    if RNGCCT(iT) == RNGCCTPL(iTPL)
        figure(FIGCCT);
        subplot(FIGCCTPLNR,FIGCCTPLNC,iTPL);
        [x,y] = planckXY(RNGCCT(iT));
        [S,R,G,B,tr,tg,tb] = RGB.getPSD(x,y);
        plot(R.npsd.X,(tr/S.npsd.Ymax)*R.npsd.Y,PLTXLCS{1});
        hold on;
        plot(G.npsd.X,(tg/S.npsd.Ymax)*G.npsd.Y,PLTXLCS{2});
        plot(B.npsd.X,(tb/S.npsd.Ymax)*B.npsd.Y,PLTXLCS{3});
        axis([FIGLDAMIN FIGLDAMAX 0 1]);
        xlabel('Wavelength (nm)');
        ylabel('Normalized PSD');
        title(sprintf('CCT = %dK',RNGCCTPL(iTPL)));
        iTPL = iTPL+1;
    end
    
    FIGBER(iT) = figure('Name',sprintf('BER vs SNR, CCT = %dK, Illumination = %dlx',RNGCCT(iT),lkIl),'NumberTitle',FIGTITLE);
    for iOf = 1:LENOFDMTYPES
        ofdmType = lower(RNGOFDMTYPES{iOf});
        switch lower(ofdmType)
            case 'acoofdm'
                STRTITLE = sprintf('ACO-OFDM, CCT = %dK, Illumination = %dlx',RNGCCT(iT),lkIl);
                plLC = PLACOLC; plLS = PLACOLS; plMK = PLACOMK; plLGD = 'ACO';
                subplot(FIGBERPLNR,FIGBERPLNC,iOf);
            case 'dcoofdm'
                STRTITLE = sprintf('DCO-OFDM, CCT = %dK, Illumination = %dlx',RNGCCT(iT),lkIl);
                plLC = PLDCOLC; plLS = PLDCOLS; plMK = PLDCOMK; plLGD = 'DCO';
                subplot(FIGBERPLNR,FIGBERPLNC,iOf);
        end
        [Xp,Yp] = getCleanPoints(RNGSNRDB,log10(sum(BER(:,:,iT,iOf),1)/NTX),PLOTDMIN);
        semilogy(Xp,Yp,[plLC plLS plMK]);
        hold on;
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        for iTx = 1:NTX
            [Xp,Yp] = getCleanPoints(RNGSNRDB,log10(BER(iTx,:,iT,iOf)),PLOTDMIN);
            semilogy(Xp,power(10,Yp),[PLTXLCS{iTx} PLTXLSS{iTx} PLTXMKS{iTx}]);
        end
        grid on;
        FIGBERLGD(:,iOf) = {plLGD;[plLGD ':Red'];[plLGD ':Green'];[plLGD ':Blue']};
        legend(gca,FIGBERLGD(:,iOf));
        xlabel('SNR^{tx}_{avg}');
        ylabel('BER');
        title(STRTITLE);
        
        LOOPCOUNT = LOOPCOUNT+1;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        waitbar(PROGRESS,h,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));
    end
end

% Figure filter responses
FIGFILT = figure('Name',sprintf('Filter Transmission'),'NumberTitle',FIGTITLE);
plot(Rrx.sensor.filter.X,Rrx.sensor.filter.Y,'r');
hold on;
plot(Grx.sensor.filter.X,Grx.sensor.filter.Y,'g');
plot(Brx.sensor.filter.X,Brx.sensor.filter.Y,'b');
axis([FIGLDAMIN FIGLDAMAX 0 1]);
xlabel('Wavelength (nm)');
ylabel('Transmission');
title('Filter Transmission');
LOOPCOUNT = LOOPCOUNT+1;
PROGRESS = LOOPCOUNT/TOTALLOOPS;
waitbar(PROGRESS,h,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));

% Figure receiver responsivisities
FIGRESP = figure('Name',sprintf('Receiver Responsivity'),'NumberTitle',FIGTITLE);
plot(Rrx.sensor.responsivity.X,Rrx.sensor.responsivity.Y,'r--');
hold on;
plot(Grx.sensor.responsivity.X,Grx.sensor.responsivity.Y,'g:');
plot(Brx.sensor.responsivity.X,Brx.sensor.responsivity.Y,'b-.');
axis([FIGLDAMIN FIGLDAMAX 0 1]);
xlabel('Wavelength (nm)');
ylabel('Responsivity (A.W^{-1})');
title('Receiver responsivity');
LOOPCOUNT = LOOPCOUNT+1;
PROGRESS = LOOPCOUNT/TOTALLOOPS;
waitbar(PROGRESS,h,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));

% Plot illuminance
% FIGILL = figure('Name',sprintf('Illuminance'),'NumberTitle',FIGTITLE);
room.drawIlluminance;
FIGILL = gcf;
set(FIGILL,'Name',sprintf('Illuminance'),'NumberTitle',FIGTITLE);
LOOPCOUNT = LOOPCOUNT+1;
PROGRESS = LOOPCOUNT/TOTALLOOPS;
waitbar(PROGRESS,h,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));

%% Save Figures
LOOPCOUNT = 0;
TOTALLOOPS = 4+LENCCT;
if fSAVEALL
    f = figure(FIGILL);
    fname = [ctDirRes STRPREFIX 'Illumination' CHARIDXARCHIVE];
    saveas(f,[fname '.png'],'png');
    saveas(f,[fname '.fig'],'fig');
    saveas(f,[fname '.eps'],'epsc');
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,h,sprintf('Saving Results: %0.2f%% done...',PROGRESS*100));
    
    f = figure(FIGFILT);
    fname = [ctDirRes STRPREFIX 'FiltTrans' CHARIDXARCHIVE];
    saveas(f,[fname '.png'],'png');
    saveas(f,[fname '.fig'],'fig');
    saveas(f,[fname '.eps'],'epsc');
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,h,sprintf('Saving Results: %0.2f%% done...',PROGRESS*100));
    
    f = figure(FIGRESP);
    fname = [ctDirRes STRPREFIX 'RecvResp' CHARIDXARCHIVE];
    saveas(f,[fname '.png'],'png');
    saveas(f,[fname '.fig'],'fig');
    saveas(f,[fname '.eps'],'epsc');
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,h,sprintf('Saving Results: %0.2f%% done...',PROGRESS*100));
    
    f = figure(FIGCCT);
    fname = [ctDirRes STRPREFIX 'SPDs' CHARIDXARCHIVE];
    saveas(f,[fname '.png'],'png');
    saveas(f,[fname '.fig'],'fig');
    saveas(f,[fname '.eps'],'epsc');
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,h,sprintf('Saving Results: %0.2f%% done...',PROGRESS*100));
    
    for iT = 1:LENCCT
        STRCCT = sprintf('%dK',RNGCCT(iT));
        f = figure(FIGBER(iT));
        fname = [ctDirRes STRPREFIX 'BERvsSNR_' STRCCT CHARIDXARCHIVE];
        saveas(f,[fname '.png'],'png');
        saveas(f,[fname '.fig'],'fig');
        saveas(f,[fname '.eps'],'epsc');
        LOOPCOUNT = LOOPCOUNT+1;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        waitbar(PROGRESS,h,sprintf('Saving Results: %0.2f%% done...',PROGRESS*100));
    end
end
delete(h);
%% restore defaults
set(0,'DefaultLineMarkerSize',dlinems);
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
% set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);

end