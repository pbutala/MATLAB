% scrOSM_ofdm
close all;
clearvars;
clc;

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fARCHIVE = true;
CHAROVERWRITE = '~';
STRPREFIX = '8_SIS_IMG_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
STRPREFIXNI = '8_SIS_NI_';
CHARIDXARCHIVENI = '';

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\11b. SISOFDM all\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/11b. SISOFDM all/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\11b. SISOFDM all\\';
    otherwise
        error('Station not defined');
end
ctFileVars = [ctDirRes STRPREFIX 'scrSISOFDM_IMG' CHARIDXARCHIVE '.mat'];      % file to store workspace
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
PLOTDMIN = 4;
SNRD = zeros(lenOfstSD,lenM);
SNRA = zeros(lenOfstSD,lenM);
% for iOfst = 1:lenOfstSD
 for iOfst = lenOfstSD
    STROFST = sprintf('OfstSD_%0.1f_%0.1f',rngOfstSDDco(iOfst),rngOfstSDAco(iOfst));
%     STROFST = sprintf('OfstSD_%0.1f',rngOfstSDDco(iOfst));
    for iM = 1:lenM
        STRM = sprintf('_adM_%d_%d',rngMaco(iM),rngMdco(iM));
        clLgdall = {};
        clLgdSplitall = {};
        
        figBER(iOfst,iM) = figure;
        set(gca,'YScale','log');
        hold all;
        fnameNI = [ctDirResNI STRPREFIXNI STROFST STRM '_BER vs SNR' CHARIDXARCHIVENI];
        uiopen([fnameNI '.fig'],1);
        figBERall(iOfst,iM) = gcf;
        [~,~,~,text_strings] = legend;
        for l=1:numel(text_strings)
            clLgdall{end+1} = text_strings{l};
        end
        figBERSplit(iOfst,iM) = figure;
        set(gca,'YScale','log');
        hold all;
        fnameNI = [ctDirResNI STRPREFIXNI STROFST STRM '_BER vs SNR Split' CHARIDXARCHIVENI];
        uiopen([fnameNI '.fig'],1);
        figBERSplitall(iOfst,iM) = gcf;
        [~,~,~,text_strings] = legend;
        for l=1:numel(text_strings)
            clLgdSplitall{end+1} = text_strings{l};
        end
        
        clLgd = {};
        clLgdSplit = {};
        iLC = 0;
        iLS = 0;
        iMK = 0;
        for iOfdm = 1:lenOfdmType
            ofdmType = rngOfdmType{iOfdm};
            iSnrTh = find(bit_err(:,iM,iOfdm,iOfst) <= BERTH,1,'first');
            switch lower(ofdmType)
                case 'acoofdm'
                    STROFDM = 'ACO';
                    STRM = sprintf('M:%d',rngMaco(iM));
                    if ~isempty(iSnrTh)
                        SNRA(iOfst,iM) = rngSNRdb(iSnrTh);
                    else
                        SNRA(iOfst,iM) = NaN;
                    end
                case {'dcoofdm', 'dmt'}
                    STROFDM = 'DCO';
                    STRM = sprintf('M:%d',rngMdco(iM));
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
            clLgd{end+1} = [STRIM ' ' STROFDM ' ' STRM];
            clLgdall{end+1} = clLgd{end};
            legend(gca,clLgd,'Location','NorthEast');
            xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
            ylabel('BER');
            tStr = sprintf('BER vs SNR, Offset: DCO=%0.1f SD ACO=%0.1f SD\nN_{sc}:%d, DCO:%d/ACO:%d bits/sym',rngOfstSDDco(iOfst),rngOfstSDAco(iOfst),Nsc,BPS(iM,1),BPS(iM,2));
            title(tStr);
            grid on;
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 0.9*BERTH 1]);
            
            set(0,'CurrentFigure',figBERall(iOfst,iM));
            iLS = rem(2,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            legend(gca,clLgdall,'Location','NorthEast');
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 0.9*BERTH 1]);
            
            set(0,'CurrentFigure',figBERSplit(iOfst,iM));
            iLS = rem(3,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_sym(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err_sym(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            clLgdSplit{end+1} = [STRIM ' ' STROFDM ' ' STRM '- SYM'];
            clLgdSplitall{end+1} = clLgdSplit{end};
            iLS = rem(4,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_led(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err_led(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            clLgdSplit{end+1} = [STRIM ' ' STROFDM ' ' STRM '- TX'];
            clLgdSplitall{end+1} = clLgdSplit{end};
            
            legend(gca,clLgdSplit,'Location','NorthEast');
            xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
            ylabel('BER (Individual)');
            tStr = sprintf('BER (Split) vs SNR, Offset: DCO=%0.1f SD ACO=%0.1f SD\nN_{sc}:%d, DCO:%d/ACO:%d bits/sym',rngOfstSDDco(iOfst),rngOfstSDAco(iOfst),Nsc,BPS(iM,1),BPS(iM,2));
            title(tStr);
            grid on;
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 0.9*BERTH 1]);
            
            set(0,'CurrentFigure',figBERSplitall(iOfst,iM));
            iLS = rem(5,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_sym(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err_sym(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            iLS = rem(6,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_led(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err_led(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            legend(gca,clLgdSplitall,'Location','NorthEast');
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 0.9*BERTH 1]);
        end % Ofdm
    end % M
end % Ofst

% if lenOfstSD > 1
%     clLgdofstAll = {};
%     fnameNIofst = [ctDirResNI STRPREFIXNI 'SNR vs Offset' CHARIDXARCHIVENI];
%     uiopen([fnameNIofst '.fig'],1);
%     figOFSTNI = gcf;
%     [~,~,~,text_strings] = legend;
%     for l=1:numel(text_strings)
%         clLgdofstAll{end+1} = text_strings{l};
%     end
%         
%     figOFST = figure;
%     clLgdOfst = {};
%     hold all;
%     iMK = 0;
%     for iM = 1:lenM
%         iLC = rem(1,lenLC)+1;
%         iLS = rem(1,lenLS)+1;
%         iMK = rem(iMK+1,lenMK)+1;
%         plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%         plot(rngOfstSDDco(1:end-1),SNRD(1:end-1,iM)-150,plStyle);
%         axis([rngSNRdb(1)-150 rngSNRdb(end)-150 0.9*BERTH 1]);
%         STRM = sprintf('M:%d',rngMdco(iM));
%         clLgdOfst{end+1} = [STRIM ' ' 'DCO ' STRM];
%         
%         iLC = rem(2,lenLC)+1;
%         iMK = rem(iMK+1,lenMK)+1;
%         plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
% %         plot(rngOfstSDDco,SNRA(:,iM),plStyle);
%         plot(rngOfstSDAco(1:end-1),SNRA(1:end-1,iM)-150,plStyle);
% 
%         STRM = sprintf('M:%d',rngMaco(iM));
%         clLgdOfst{end+1} = [STRIM ' ' 'ACO ' STRM];
%     end
%     xlabel('Offset (Std Dev)');
%     yStr = sprintf('{SNR_{avg}^{tx}} - 150 (dB)');
%     ylabel(yStr);
%     tStr = sprintf('SNR vs Offset, Target BER = 10^{%d}',log10(BERTH));
%     title(tStr);
%     grid on;
% %     axis([rngOfstSDDco(1) rngOfstSDDco(end) rngSNRdb(1)-150 rngSNRdb(end)-150]);
%     legend(gca,clLgdOfst,'Location','NorthEast');
%     
%     set(0,'CurrentFigure',figOFSTNI);
%     hold all;
%     for iM = 1:lenM
%         iLC = rem(1,lenLC)+1;
%         iLS = rem(1,lenLS)+1;
%         iMK = rem(iMK+1,lenMK)+1;
%         plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%         plot(rngOfstSDDco(1:end-1),SNRD(1:end-1,iM)-150,plStyle);
%         STRM = sprintf('M:%d',rngMdco(iM));
%         clLgdofstAll{end+1} = [STRIM ' ' 'DCO ' STRM];
%         
%         iLC = rem(2,lenLC)+1;
%         iMK = rem(iMK+1,lenMK)+1;
%         plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
% %         plot(rngOfstSDDco,SNRA(:,iM),plStyle);
%         plot(rngOfstSDAco(1:end-1),SNRA(1:end-1,iM)-150,plStyle);
% 
%         STRM = sprintf('M:%d',rngMaco(iM));
%         clLgdofstAll{end+1} = [STRIM ' ' 'ACO ' STRM];
%     end
%     xlabel('Offset (Std Dev)');
%     yStr = sprintf('{SNR_{avg}^{tx}} - 150 (dB)');
%     ylabel(yStr);
%     tStr = sprintf('SNR vs Offset, Target BER = 10^{%d}',log10(BERTH));
%     title(tStr);
%     grid on;
%     axis([rngOfstSDDco(1) rngOfstSDDco(end-1) rngSNRdb(1)-150 rngSNRdb(end)-150]);
%     legend(gca,clLgdofstAll,'Location','NorthEast');
% end

% % draw setup
% figSetup = figure;
% room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
% rotate3d on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fSAVEALL
%     copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
%     copyfile(ctFileCodeSrcNI,ctFileCodeDestNI); % save script NI
% %     save(ctFileVars);                       % save workspace
%     fname = [ctDirRes STRPREFIX 'Setup' CHARIDXARCHIVE];
%     saveas(figSetup,[fname '.png'],'png');
%     saveas(figSetup,[fname '.fig'],'fig');
%     saveas(figSetup,[fname '.eps'],'epsc');
%     for iOfst = 1:lenOfstSD
    for iOfst = lenOfstSD
        STROFST = sprintf('OfstSD_%0.1f_%0.1f',rngOfstSDDco(iOfst),rngOfstSDAco(iOfst));
        for iM = 1:lenM
            STRM = sprintf('_adM_%d_%d',rngMaco(iM),rngMdco(iM));
            f = figure(figBER(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRM '_BER vs SNR' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
            f = figure(figBERall(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRM '_ALL BER vs SNR' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
            
            f = figure(figBERSplit(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRM '_BER vs SNR Split' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
            f = figure(figBERSplitall(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRM '_ALL BER vs SNR Split' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
        end
    end
%     if lenOfstSD > 1
%         fname = [ctDirRes STRPREFIX 'SNR vs Offset' CHARIDXARCHIVE];
%         saveas(figOFST,[fname '.png'],'png');
%         saveas(figOFST,[fname '.fig'],'fig');
%         saveas(figOFST,[fname '.eps'],'epsc');
%         fname = [ctDirRes STRPREFIX 'ALL SNR vs Offset' CHARIDXARCHIVE];
%         saveas(figOFSTNI,[fname '.png'],'png');
%         saveas(figOFSTNI,[fname '.fig'],'fig');
%         saveas(figOFSTNI,[fname '.eps'],'epsc');
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





































