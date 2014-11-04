% scrCSKOFDMPL
if ~exist('ctFileVars','var')
    error('Data file not specified');
end
if (exist('hWB','var')&& ishandle(hWB))
    delete(hWB);
end

close all;
clearvars -except ctFileVars;
clc;
load(ctFileVars);                                                             % LOAD DATA
ctFilePlot = [ctDirData STRPREFIX 'plCSKOFDM' CHARIDXARCHIVE '.mat'];   % Plot file name

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');
dfigppm = get(0,'DefaultFigurePaperPositionMode');
set(0,'DefaultFigurePaperPositionMode','Manual');
dfigpu = get(0,'DefaultFigurePaperUnits');
set(0,'DefaultFigurePaperUnits','Inches');
dfigpp = get(0,'DefaultFigurePaperPosition');
set(0,'DefaultFigurePaperPosition',[0 0 8 6]);
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',8);
FIGTITLE = 'Off';

try
    % Wait Bar to show progress
    if fSHOWPGBAR
        hWB = waitbar(0,'Plotting Results: 0.00% done','Name',WBTITLE,'Visible','Off');
        set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
    end
    % PLOT Configs
    PLOTDMIN = 5;
    RNGSNROFST = RNGSNRDB - SNROFST;
    FIGBERXMIN = RNGSNRMIN-SNROFST; FIGBERXMAX = RNGSNRMAX-SNROFST;
    SNRatTBER = nan(LENMOD,LENMODNSC,LENOFDMOFST,LENOFDMTYPES);
    
    PLACOLC = 'm'; PLACOLS = '-'; PLACOMK = 'o';
    PLDCOLC = 'c'; PLDCOLS = '-'; PLDCOMK = 'd';
    PLTXLCS = {'r';'g';'b'}; PLTXLSS = {'--';'-.';':'}; PLTXMKS = {'>';'s';'*'};
    PLNSCMKS = {'x';'h';'^';'+';'v';'*';'<';'p'};
    MKTYP = {'o','+','^','s'};
    MKCLR = {'g','y','b','r'};
    O_LC = {'m';'c'}; O_LS = {'--';':'}; O_MK = {'+','^'};
    O_MKsp_O = {'x','o'}; O_MKsp_C = {'*','s'};
    % MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};
    % Figure BER vs SNR config
    FIGBERYMIN = 0.9*BERTH; FIGBERYMAX = 1;
    
    LOOPCOUNT = 0;
    TOTALLOOPS = 2 + 2*LENMOD*LENOFDMOFST*LENMODNSC + LENMOD*LENMODNSC;
    
    % ********************** PLOT AND SAVE SPDs responsivities *********************
    FIGSPD = figure('Name','Spectral Spreads','NumberTitle',FIGTITLE);
    hold on;
    sLGD = {}; hLGD = [];
    for i=1:RGB.NCLR
        p = RGB.PSDs{i};
        LSTL = [PLTXLCS{i} PLTXLSS{i}]; MSTL = [PLTXLCS{i} PLTXMKS{i}]; PSTL = [PLTXLCS{i} PLTXLSS{i} PLTXMKS{i}];
        Yp = p.npsd.Ymax; Xp = p.npsd.X(p.npsd.Y == Yp);
        plot(p.npsd.X, p.npsd.Y, LSTL);
        plot(Xp, Yp, MSTL);
        hLGD(end+1) = plot(nan,nan,PSTL);
        sLGD{end+1} = ['SPD: ' p.name];
    end
    V = getEyeSens(LAMBDAMIN,LAMBDAMAX,LAMBDADELTA,1978);
    [Xp,Yp] = getCleanPoints(lambdas,V/max(V), 50);
    plot(lambdas, V/max(V), 'k-');
    plot(Xp, Yp, ['k' PLNSCMKS{1}]);
    hLGD(end+1) = plot(nan,nan,['k-' PLNSCMKS{1}]);
    sLGD{end+1} = 'Eye Sensitivity';
    
    legend(gca,hLGD,sLGD);
    grid on;
    
    if fSAVEALL
        fname = [ctDirData STRPREFIX 'SpectralSpreads' CHARIDXARCHIVE];
        saveas(FIGSPD,[fname '.png'],'png');
        saveas(FIGSPD,[fname '.fig'],'fig');
        saveas(FIGSPD,[fname '.eps'],'epsc');
    end
    if fCLOSEALL
        close(FIGSPD);
    end
    % Update Wait bar
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    if fSHOWPGBAR
        waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
    end
    
    % ********************** PLOT AND SAVE COLOR GAMUT*********************
    FIGGMT = figure('Name','RGB LED xy gamut','NumberTitle',FIGTITLE);
    RGB.obs.showGamut();
    hold on;
    scatter(RGB.xyz(1,:),RGB.xyz(2,:),8,'.','y');
    [x,y,~] = RGB.obs.getCoordinates(RGB.PSDs{1});
    scatter(x,y,80,'x','r','LineWidth',2);
    [x,y,~] = RGB.obs.getCoordinates(RGB.PSDs{2});
    scatter(x,y,80,'x','g','LineWidth',2);
    [x,y,~] = RGB.obs.getCoordinates(RGB.PSDs{3});
    scatter(x,y,80,'x','b','LineWidth',2);
    
    title(sprintf('RGB LED xy gamut\nR:%dnm G:%dnm B:%dnm Res:%0.2f',RMN,GMN,BMN,RES));
    if fSAVEALL
        fname = [ctDirData STRPREFIX 'Gamut' CHARIDXARCHIVE];
        saveas(FIGGMT,[fname '.png'],'png');
        saveas(FIGGMT,[fname '.fig'],'fig');
        saveas(FIGGMT,[fname '.eps'],'epsc');
    end
    if fCLOSEALL
        close(FIGGMT);
    end
    % Update Wait bar
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    if fSHOWPGBAR
        waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
    end
    
    % ********************** PLOT AND SAVE BER vs SNR *********************
    for iM = 1:LENMOD
        for iNSC = 1:LENMODNSC                                                  % LOOP START MODNSC
            MODNSC = RNGMODNSC(iNSC);
            STRNSC = sprintf(' Nsc %d', MODNSC);
            for iDC = 1:LENOFDMOFST
                STRDCOOFST = sprintf(' D %0.2f SD', RNGOFDMOFSTDCO(iDC));
                STRACOOFST = sprintf(' A %0.2f SD', RNGOFDMOFSTACO(iDC));
                
                FIGBER(iNSC,iDC) = figure;
                set(gca,'YScale','log');
                hold on;
                sLGD = {}; hLGD = [];
                for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
                    ofdmType = lower(RNGOFDMTYPES{iOf});
                    switch lower(ofdmType)
                        case 'acoofdm'
                            STROFDM = 'ACO';
                            STRMODACO = sprintf(' A M-%d',RNGMODACO(iM));
                        case 'dcoofdm'
                            STROFDM = 'DCO';
                            STRMODDCO = sprintf(' D M-%d',RNGMODDCO(iM));
                    end
                    [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iM,iNSC,iDC,iOf),log10(BER(:,iM,iNSC,iDC,iOf)),PLOTDMIN);  % Get points well spaced out
                    semilogy(Xp,power(10,Yp),[O_LC{iOf} O_MK{iOf}]);                          % Semilog AVG BER vs SNR Marker
                    semilogy(RNGSNROFST(:,iM,iNSC,iDC,iOf),BER(:,iM,iNSC,iDC,iOf),[O_LC{iOf} O_LS{iOf}]);  % Semilog AVG BER vs SNR
                    hLGD(end+1) = semilogy(nan,nan,[O_LC{iOf} O_LS{iOf} O_MK{iOf}]);
                    sLGD{end+1} = ['CSK-' STROFDM];
                    
                    %** START: find SNR for target BER **
                    iTH = find(BER(:,iM,iNSC,iDC,iOf) <= BERTH,1,'first');
                    if ~isempty(iTH)
                        SNRatTBER(iM,iNSC,iDC,iOf) = RNGSNROFST(iTH,iM,iNSC,iDC,iOf);
                    end
                    %** END: find SNR for target BER **
                end
                axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
                grid on;
                xlabel(sprintf('SNR^{tx}_{avg}(dB)'));
                ylabel('BER');
                strTitle = sprintf(['DCO: %3d BPS | %3d-QAM | %0.2f SD\n'...
                                    'ACO: %3d BPS | %3d-QAM | %0.2f SD\n'...
                                    '%2d-CSK: %3d BPS | %3d N_{sc}'],...
                                    O_BPS(iM,iNSC,1), RNGMODDCO(iM), RNGOFDMOFSTDCO(iDC),...
                                    O_BPS(iM,iNSC,2), RNGMODACO(iM), RNGOFDMOFSTACO(iDC),...
                                    C_M, C_BPS(iNSC), MODNSC);
                title(strTitle);
                legend(gca,hLGD,sLGD);
                
                f = gcf;
                if fSAVEALL
                    fname = [ctDirData STRPREFIX 'BERvsSNR' STRMODDCO STRMODACO STRNSC STRDCOOFST STRACOOFST CHARIDXARCHIVE];
                    saveas(f,[fname '.png'],'png');
                    saveas(f,[fname '.fig'],'fig');
                    saveas(f,[fname '.eps'],'epsc');
                end
                if fCLOSEALL
                    close(f);
                end
                % Update Wait bar
                LOOPCOUNT = LOOPCOUNT+1;
                PROGRESS = LOOPCOUNT/TOTALLOOPS;
                if fSHOWPGBAR
                    waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
                end
            end
        end
    end
    
    % ******************* PLOT AND SAVE SPLIT BER vs SNR ******************
    for iM = 1:LENMOD
        for iNSC = 1:LENMODNSC                                                  % LOOP START MODNSC
            MODNSC = RNGMODNSC(iNSC);
            STRNSC = sprintf(' Nsc %d', MODNSC);
            for iDC = 1:LENOFDMOFST
                STRDCOOFST = sprintf(' D %0.2f SD', RNGOFDMOFSTDCO(iDC));
                STRACOOFST = sprintf(' A %0.2f SD', RNGOFDMOFSTACO(iDC));
                
                FIGBERSPT(iNSC,iDC) = figure;
                set(gca,'YScale','log');
                hold on;
                sLGD = {}; hLGD = [];
                for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
                    ofdmType = lower(RNGOFDMTYPES{iOf});
                    switch lower(ofdmType)
                        case 'acoofdm'
                            STROFDM = 'ACO';
                            STRMODACO = sprintf(' A M-%d',RNGMODACO(iM));
                        case 'dcoofdm'
                            STROFDM = 'DCO';
                            STRMODCO = sprintf(' D M-%d',RNGMODDCO(iM));
                    end
                    [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iM,iNSC,iDC,iOf),log10(C_BER(:,iM,iNSC,iDC,iOf)),PLOTDMIN);  % Get points well spaced out
                    semilogy(Xp,power(10,Yp),[O_LC{iOf} O_MKsp_C{iOf}]);                          % Semilog AVG BER vs SNR Marker
                    semilogy(RNGSNROFST(:,iM,iNSC,iDC,iOf),C_BER(:,iM,iNSC,iDC,iOf),[O_LC{iOf} O_LS{iOf}]);  % Semilog AVG BER vs SNR
                    hLGD(end+1) = semilogy(nan,nan,[O_LC{iOf} O_LS{iOf} O_MKsp_C{iOf}]);
                    sLGD{end+1} = ['CSK (' STROFDM ')'];
                    
                    [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iM,iNSC,iDC,iOf),log10(O_BER(:,iM,iNSC,iDC,iOf)),PLOTDMIN);  % Get points well spaced out
                    semilogy(Xp,power(10,Yp),[O_LC{iOf} O_MKsp_O{iOf}]);                          % Semilog AVG BER vs SNR Marker
                    semilogy(RNGSNROFST(:,iM,iNSC,iDC,iOf),O_BER(:,iM,iNSC,iDC,iOf),[O_LC{iOf} O_LS{iOf}]);  % Semilog AVG BER vs SNR
                    hLGD(end+1) = semilogy(nan,nan,[O_LC{iOf} O_LS{iOf} O_MKsp_O{iOf}]);
                    sLGD{end+1} = ['OFDM (' STROFDM ')'];
                end
                axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
                grid on;
                xlabel(sprintf('SNR^{tx}_{avg}(dB)'));
                ylabel('BER');
                strTitle = sprintf(['DCO: %3d BPS | %3d-QAM | %0.2f SD\n'...
                                    'ACO: %3d BPS | %3d-QAM | %0.2f SD\n'...
                                    '%2d-CSK: %3d BPS | %3d N_{sc}'],...
                                    O_BPS(iM,iNSC,1), RNGMODDCO(iM), RNGOFDMOFSTDCO(iDC),...
                                    O_BPS(iM,iNSC,2), RNGMODACO(iM), RNGOFDMOFSTACO(iDC),...
                                    C_M, C_BPS(iNSC), MODNSC);
                title(strTitle);
                legend(gca,hLGD,sLGD);
                
                f = gcf;
                if fSAVEALL
                    fname = [ctDirData STRPREFIX 'BER_SplitvsSNR' STRMODDCO STRMODACO STRNSC STRDCOOFST STRACOOFST CHARIDXARCHIVE];
                    saveas(f,[fname '.png'],'png');
                    saveas(f,[fname '.fig'],'fig');
                    saveas(f,[fname '.eps'],'epsc');
                end
                if fCLOSEALL
                    close(f);
                end
                % Update Wait bar
                LOOPCOUNT = LOOPCOUNT+1;
                PROGRESS = LOOPCOUNT/TOTALLOOPS;
                if fSHOWPGBAR
                    waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
                end
            end
        end
    end
    
    % ********************** PLOT AND SAVE SNR vs DC OFSTs *********************
    if LENOFDMOFST > 1
        for iM = 1:LENMOD
            for iNSC = 1:LENMODNSC                                                  % LOOP START MODNSC
                MODNSC = RNGMODNSC(iNSC);
                STRNSC = sprintf(' Nsc %d', MODNSC);
                
                FIGSNRvsDC(iNSC,iDC) = figure;
                hold on;
                sLGD = {}; hLGD = [];
                for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
                    ofdmType = lower(RNGOFDMTYPES{iOf});
                    switch lower(ofdmType)
                        case 'acoofdm'
                            STROFDM = 'ACO';
                            RNGDC = RNGOFDMOFSTACO;
                            STRMODACO = sprintf(' A M-%d',RNGMODACO(iM));
                        case 'dcoofdm'
                            STROFDM = 'DCO';
                            RNGDC = RNGOFDMOFSTDCO;
                            STRMODDCO = sprintf(' D M-%d',RNGMODDCO(iM));
                    end
                    if min(RNGDC(:)) ~= max(RNGDC(:));
                        [Xp,Yp] = getCleanPoints(RNGDC(1:end-LENOFSTIGNR),SNRatTBER(iM,iNSC,1:end-LENOFSTIGNR,iOf),DOFST);  % Get points well spaced out
                        plot(Xp,Yp,[O_LC{iOf} O_MK{iOf}]);                          % Semilog AVG BER vs SNR Marker
                        plot(RNGDC(1:end-LENOFSTIGNR),squeeze(SNRatTBER(iM,iNSC,1:end-LENOFSTIGNR,iOf)),[O_LC{iOf} O_LS{iOf}]);  % Semilog AVG BER vs SNR
                        hLGD(end+1) = plot(nan,nan,[O_LC{iOf} O_LS{iOf} O_MK{iOf}]);
                        sLGD{end+1} = ['CSK-' STROFDM];
                    end
                end
                axis([min(RNGDC) max(RNGDC) FIGBERXMIN FIGBERXMAX]);
                grid on;
                xlabel(sprintf('DC Offset (SD)'));
                ylabel('SNR^{tx}_{avg}(dB)');
                strTitle = sprintf(['DCO: %3d BPS | %3d-QAM | '...
                    'ACO: %3d BPS | %3d-QAM\n'...
                    '%2d-CSK: %3d BPS | %3d N_{sc}'],...
                    O_BPS(iM,iNSC,1), RNGMODDCO(iM),...
                    O_BPS(iM,iNSC,2), RNGMODACO(iM),...
                    C_M, C_BPS(iNSC),MODNSC);
                title(strTitle);
                legend(gca,hLGD,sLGD);
                
                f = gcf;
                if fSAVEALL
                    fname = [ctDirData STRPREFIX 'SNRvsDCOFST' STRMODDCO STRMODACO STRNSC CHARIDXARCHIVE];
                    saveas(f,[fname '.png'],'png');
                    saveas(f,[fname '.fig'],'fig');
                    saveas(f,[fname '.eps'],'epsc');
                end
                if fCLOSEALL
                    close(f);
                end
                % Update Wait bar
                LOOPCOUNT = LOOPCOUNT+1;
                PROGRESS = LOOPCOUNT/TOTALLOOPS;
                if fSHOWPGBAR
                    waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
                end
            end
        end
    else
        % Update Wait bar
        LOOPCOUNT = LOOPCOUNT+LENMOD*LENMODNSC;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
        end
    end
    % ********************** ************************ *********************
    save(ctFilePlot);                       % save workspace
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
catch ex
    save(ctFilePlot);                       % save workspace
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    %% restore defaults
    set(0,'DefaultLineMarkerSize',dlinems);
    set(0,'DefaultLineLineWidth',dlinelw);
    set(0,'DefaultAxesFontName',daxesfontname);
    set(0,'DefaultAxesFontSize',daxesfontsize);
    set(0,'DefaultFigureVisible',dfigvis);
    set(0,'DefaultFigurePaperPosition',dfigpp);
    set(0,'DefaultFigurePaperUnits',dfigpu);
    set(0,'DefaultFigurePaperPositionMode',dfigppm);
    setpref('Internet','E_mail','pbutala@bu.edu');
    setpref('Internet','SMTP_Server','smtp.bu.edu');
    STREMAIL = ['Simulation ' STRPREFIX ' done with errors.'];
    sendmail('pankil.butala@gmail.com',STREMAIL);
    rethrow(ex);
end
%% restore defaults
set(0,'DefaultLineMarkerSize',dlinems);
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);
setpref('Internet','E_mail','pbutala@bu.edu');
setpref('Internet','SMTP_Server','smtp.bu.edu');
STREMAIL = ['Simulation ' STRPREFIX ' done.'];
sendmail('pankil.butala@gmail.com',STREMAIL);
% end