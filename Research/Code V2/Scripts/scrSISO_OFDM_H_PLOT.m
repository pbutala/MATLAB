% scrOSM_ofdm_NI
close all;
clearvars;
clc;

fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fARCHIVE = true;

CHAROVERWRITE = '~';
STRPREFIX = '3_SISO_H_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
STRPREFIXIMG = '3_SIS_IMG_';
CHARIDXARCHIVEIMG = '';

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\11. SISO OFDM H\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/11. SISO OFDM H/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\11. SISO OFDM H\\';
    otherwise
        error('Station not defined');
end
ctFileVars = [ctDirRes STRPREFIX 'scrSISO_OFDM_H' CHARIDXARCHIVE '.mat'];      % file to store workspace
load(ctFileVars);

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
set(0,'DefaultFigurePaperPosition',[0 0 8 6]);
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTDMIN = 4;
SNRD = zeros(lenOfstSD,lenM);
SNRA = zeros(lenOfstSD,lenM);
for iOfst = 1:lenOfstSD
    STROFSTS = sprintf('OfstSD_%0.2f_%0.2f',rngOfstSDDcoS(iOfst),rngOfstSDAcoS(iOfst));
    STROFSTM = sprintf('OfstSD_%0.2f_%0.2f',rngOfstSDDcoM(iOfst),rngOfstSDAcoM(iOfst));
    for iM = 1:lenM
        STRMS = sprintf('_adM_%d_%d',rngMacoS(iM),rngMdcoS(iM));
        STRMM = sprintf('_adM_%d_%d',rngMacoM(iM),rngMdcoM(iM));
        clLgdall = {};
        figBER(iOfst,iM) = figure;
        set(gca,'YScale','log');
        hold all;
%         figBERSplit(iOfst,iM) = figure;
%         set(gca,'YScale','log');
%         hold all;
        fnameIMG = [ctDirResIMG STRPREFIXIMG STROFSTM STRMM '_ALL BER vs SNR' CHARIDXARCHIVEIMG];
        uiopen([fnameIMG '.fig'],1);
        figBERall(iOfst,iM) = gcf;
        [~,~,~,text_strings] = legend;
        for l=1:numel(text_strings)
            clLgdall{end+1} = text_strings{l};
        end
        
        clLgd = {};
%         clLgdSplit = {};
        iLC = 0;
        iLS = 0;
        iMK = 0;
        for iOfdm = 1:lenOfdmType
            ofdmType = rngOfdmType{iOfdm};
            iSnrTh = find(bit_err(:,iM,iOfdm,iOfst) <= BERTH,1,'first');
            switch lower(ofdmType)
                case 'acoofdm'
                    STROFDM = 'ACO';
                    STRMS = sprintf('M:%d',rngMacoS(iM));
                    if ~isempty(iSnrTh)
                        SNRA(iOfst,iM) = rngSNRdb(iSnrTh);
                    else
                        SNRA(iOfst,iM) = NaN;
                    end
                case {'dcoofdm', 'dmt'}
                    STROFDM = 'DCO';
                    STRMS = sprintf('M:%d',rngMdcoS(iM));
                    if ~isempty(iSnrTh)
                        SNRD(iOfst,iM) = rngSNRdb(iSnrTh);
                    else
                        SNRD(iOfst,iM) = NaN;
                    end
            end
            
            iLC = rem(iOfdm,lenLC)+1;
            iLS = rem(iM,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            
            set(0,'CurrentFigure',figBER(iOfst,iM));
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            clLgd{end+1} = [STRSS ' ' STROFDM ' ' STRMS];
            clLgdall{end+1} = clLgd{end};
            legend(gca,clLgd,'Location','NorthEast');
%             xlabel('{P_{avg}^{tx}}/{N_{0}} (dB)');
            xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
            ylabel('BER');
            tStr = sprintf('BER vs SNR, N_{sc}:%d, Offset SISO DCO:%0.2f/ACO:%0.2f SD,\nDCO:%d/ACO:%d bits/sym',Nsc,rngOfstSDDcoS(iOfst),rngOfstSDAcoS(iOfst),BPSS(iM,1),BPSS(iM,2));
            title(tStr);
            grid on;
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 BERTH/5 1]);
            
            set(0,'CurrentFigure',figBERall(iOfst,iM));
            iLS = rem(3,lenLS)+1;
            iMK = rem(iMK+3,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            legend(gca,clLgdall,'Location','NorthEast');
            tStr = sprintf('BER vs SNR, N_{sc}:%d, Offset SISO DCO:%0.2f/ACO:%0.2f SD, SIS DCO:%0.2f/ACO:%0.2f SD\nSISO DCO:%d/ACO:%d bits/sym, SIS DCO:%d/ACO:%d bits/sym',...
                            Nsc,rngOfstSDDcoS(iOfst),rngOfstSDAcoS(iOfst),rngOfstSDDcoM(iOfst),rngOfstSDAcoM(iOfst),BPSS(iM,1),BPSS(iM,2),BPSM(iM,1),BPSM(iM,2));
            title(tStr);
%             set(0,'CurrentFigure',figBERSplit(iOfst,iM));
%             iLS = rem(2,lenLS)+1;
%             iMK = rem(iMK+1,lenMK)+1;
%             plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%             [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_sym(:,iM,iOfdm,iOfst)),PLOTDMIN);
% %             semilogy(gca,rngSNRdb,bit_err_sym(:,iM,iOfdm,iOfst),plStyle);
%             semilogy(gca,Xc-150,power(10,Yc),plStyle);
%             clLgdSplit{end+1} = [STRNI ' ' STROFDM ' ' STRM '- SYM'];
%             iLS = rem(3,lenLS)+1;
%             iMK = rem(iMK+1,lenMK)+1;
%             plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%             [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_led(:,iM,iOfdm,iOfst)),PLOTDMIN);
% %             semilogy(gca,rngSNRdb,bit_err_led(:,iM,iOfdm,iOfst),plStyle);
%             semilogy(gca,Xc-150,power(10,Yc),plStyle);
%             clLgdSplit{end+1} = [STRNI ' ' STROFDM ' ' STRM '- TX'];
%             
%             legend(gca,clLgdSplit,'Location','NorthEast');
% %             xlabel('{P_{avg}^{tx}}/{N_{0}} (dB)');
%             xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
%             ylabel('BER (Individual)');
%             tStr = sprintf('BER (Split) vs SNR, Offset: %0.2f SD\nN_{sc}:%d, DCO:%d/ACO:%d bits/sym',rngOfstSD(iOfst),Nsc,BPS(iM,1),BPS(iM,2));
%             title(tStr);
%             grid on;
%             axis([rngSNRdb(1)-150 rngSNRdb(end)-150 BERTH/5 1]);
        end % Ofdm
    end % M
end % Ofst

if lenOfstSD > 1
    figOFST = figure;
    clLgdOfst = {};
    hold all;
    iMK = 0;
    for iM = 1:lenM
        iLC = rem(1,lenLC)+1;
        iLS = rem(1,lenLS)+1;
        iMK = rem(iMK+1,lenMK)+1;
        plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
        plot(rngOfstSDDcoS,SNRD(:,iM),plStyle);
        STRMS = sprintf('M:%d',rngMdcoS(iM));
        clLgdOfst{end+1} = [STRSS ' ' 'DCO ' STRMS];
        
        iLC = rem(2,lenLC)+1;
        iMK = rem(iMK+1,lenMK)+1;
        plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
        plot(rngOfstSDDcoS,SNRA(:,iM),plStyle);
        STRMS = sprintf('M:%d',rngMacoS(iM));
        clLgdOfst{end+1} = [STRSS ' ' 'ACO ' STRMS];
    end
    xlabel('Offset (Std Dev)');
    yStr = sprintf('{SNR_{avg}^{tx}} - 150 (dB) (BER=10^{%d})',log10(BERTH));
    ylabel(yStr);
    tStr = sprintf('SNR vs Offset, Target BER = 10^{%d}',log10(BERTH));
    title(tStr);
    grid on;
    axis([rngOfstSDDcoS(1) rngOfstSDDcoS(end) rngSNRdb(1)-150 rngSNRdb(end)-150]);
    legend(gca,clLgdOfst,'Location','NorthEast');
end

% % draw setup
% figSetup = figure;
% room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
% rotate3d on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fSAVEALL
%     copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
% %     save(ctFileVars);                       % save workspace
%     fname = [ctDirRes STRPREFIX 'Setup' CHARIDXARCHIVE];
%     saveas(figSetup,[fname '.png'],'png');
%     saveas(figSetup,[fname '.fig'],'fig');
%     saveas(figSetup,[fname '.eps'],'epsc');
%     for iOfst = 1:lenOfstSD
%         STROFST = sprintf('OfstSD_%0.2f',rngOfstSDDcoS(iOfst));
%         for iM = 1:lenM
%             STRMS = sprintf('_adM_%d_%d_%d_%d',rngMacoS(iM),rngMdcoS(iM),rngMacoM(iM),rngMdcoM(iM));
%             f = figure(figBER(iOfst,iM));
%             fname = [ctDirRes STRPREFIX STROFST STRMS '_BER vs SNR' CHARIDXARCHIVE];
%             saveas(f,[fname '.png'],'png');
%             saveas(f,[fname '.fig'],'fig');
%             saveas(f,[fname '.eps'],'epsc');
%             f = figure(figBERall(iOfst,iM));
%             fname = [ctDirRes STRPREFIX STROFST STRMS '_ALL BER vs SNR' CHARIDXARCHIVE];
%             saveas(f,[fname '.png'],'png');
%             saveas(f,[fname '.fig'],'fig');
%             saveas(f,[fname '.eps'],'epsc');
% %             f = figure(figBERSplit(iOfst,iM));
% %             fname = [ctDirRes STRPREFIX STROFST STRM '_BER vs SNR Split' CHARIDXARCHIVE];
% %             saveas(f,[fname '.png'],'png');
% %             saveas(f,[fname '.fig'],'fig');
% %             saveas(f,[fname '.eps'],'epsc');
%         end
%     end
%     if lenOfstSD > 1
%         fname = [ctDirRes STRPREFIX 'SNR vs Offset' CHARIDXARCHIVE];
%         saveas(figOFST,[fname '.png'],'png');
%         saveas(figOFST,[fname '.fig'],'fig');
%         saveas(figOFST,[fname '.eps'],'epsc');
%     end
end


% restore defaults
set(0,'DefaultLineMarkerSize',dlinems);
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
% set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);

disp(' ');
disp('--DONE--');






























