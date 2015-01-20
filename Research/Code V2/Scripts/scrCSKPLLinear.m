% scrCSKPLLinear
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
    RNGSNRMINPL = 0; RNGSNRMAXPL = 25;
    PLOTDMIN = 5;
    RNGSNROFST = RNGSNRDB - SNROFST;
    FIGBERXMIN = RNGSNRMINPL - SNROFST;
    FIGBERXMAX = RNGSNRMAXPL - SNROFST;
    FIGBERYMIN = 0.9*BERTH; FIGBERYMAX = 1;
    
    % Plot optics
    PLACOLC = 'm'; PLACOLS = '-'; PLACOMK = 'o';
    PLDCOLC = 'c'; PLDCOLS = '-'; PLDCOMK = 'd';
    PLTXLCS = {'r';'g';'b'}; PLTXLSS = {'--';'-.';':'}; PLTXMKS = {'>';'s';'*'};
    PLNSCMKS = {'x';'h';'^';'+';'v';'*';'<';'p'};
    MKTYP = {'o','+','^','s'};
    MKCLR = {'g','y','b','r'};
    
    % Figure BER vs SNR config
    STRSNR = 'SNR_{avg}';
    STRLSNR = 'LSNR_{avg}';
    PLCBCLCS = {'r';'g';'b';'r';'g';'b';'r';'g';'b'};
    PLCBCLSS = {'--';'-.';':';'-.';':';'--';':';'--';'-.'};
    PLCBCMKS = {'x';'+';'*';'s';'d';'h';'<';'v';'>'};
    
    LOOPCOUNT = 0;
    TOTALLOOPS = LENCBC*2 + (LENCBC>1)*LENCBC;
    TSTART = tic;
    
    % ******************** PLOT AND SAVE CONSTELLATION ********************
    flCIE = 'CIE1931_JV_1978_2deg.csv';
    obs = cCIE(flCIE);
    for iCBC=1:LENCBC
        fCBC = RNGCBC(iCBC);
        FIGCNST = figure('Name',sprintf('CBC%d %d-CSK Constellation',fCBC,M),'NumberTitle',FIGTITLE);
        obs.showGamut();
        hold on;
        [x,y] = csk(iCBC).getSyms();
        scatter(x,y,80,'k','x','linewidth',2);
        
        title(sprintf('%d-CSK constellation for CBC%d',M,fCBC));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('%d-CSK_CBC%d_Constellation',M,fCBC) CHARIDXARCHIVE];
            saveas(FIGCNST,[fname '.png'],'png');
            saveas(FIGCNST,[fname '.fig'],'fig');
            saveas(FIGCNST,[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGCNST);
        end
        % Update Wait bar
        LOOPCOUNT = LOOPCOUNT+1;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
        end
    end
    
    % ******************** PLOT AND SAVE BER vs SNR_opt *******************
    if LENCBC > 1
        sLGD = {}; hLGD = [];
        FIGBERALL = figure('Name',sprintf('%d-CSK BER vs SNR',M),'NumberTitle',FIGTITLE);
        set(gca,'YScale','Log');
        hold on;
    end
    
    for iCBC=1:LENCBC
        fCBC = RNGCBC(iCBC);
        LSTL = [PLCBCLCS{iCBC} PLCBCLSS{iCBC}];
        MSTL = [PLCBCLCS{iCBC} PLCBCMKS{iCBC}];
        PSTL = [PLCBCLCS{iCBC} PLCBCLSS{iCBC} PLCBCMKS{iCBC}];
        FIGBER(iCBC) = figure('Name',sprintf('%d-CSK BER vs SNR CBC%d',M,fCBC),'NumberTitle',FIGTITLE);
        [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iCBC),log10(BER(:,iCBC)),PLOTDMIN);  % Get points well spaced out
        semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker
        hold on;
        semilogy(RNGSNROFST(:,iCBC),BER(:,iCBC),LSTL);          % Semilog AVG BER vs SNR
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        grid on;
        xlabel('SNR_{opt}(dB)');
        ylabel('BER');
        title(sprintf('%d-CSK BER vs SNR CBC%d',M,fCBC));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('%d-CSK_CBC%d_BERvsSNR',M,fCBC) CHARIDXARCHIVE];
            saveas(FIGBER(iCBC),[fname '.png'],'png');
            saveas(FIGBER(iCBC),[fname '.fig'],'fig');
            saveas(FIGBER(iCBC),[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGBER(iCBC));
        end
        % Update Wait bar
        LOOPCOUNT = LOOPCOUNT+1;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
        end
        
        if LENCBC > 1
            % BER vs SNR_opt
            figure(FIGBERALL);
            semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker\\
            semilogy(RNGSNROFST(:,iCBC),BER(:,iCBC),LSTL);          % Semilog AVG BER vs SNR
            hLGD(end+1) = semilogy(nan,nan,PSTL);
            sLGD{end+1} = sprintf('CBC%d',fCBC);
            
            % Update Wait bar (snr)
            LOOPCOUNT = LOOPCOUNT+1;
            PROGRESS = LOOPCOUNT/TOTALLOOPS;
            TELAPSED = toc(TSTART);
            TREM = (TELAPSED/PROGRESS)-TELAPSED;
            if fSHOWPGBAR
                waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
            end
        end
    end
    if LENCBC > 1
        % snr_opt
        figure(FIGBERALL);
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        legend(gca,hLGD,sLGD,'Location','NorthEast');
        grid on;
        xlabel('SNR_{opt}(dB)');
        ylabel('BER');
        title(sprintf('%d-CSK BER vs SNR',M));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('%d-CSK_BERvsSNR',M) CHARIDXARCHIVE];
            saveas(FIGBERALL,[fname '.png'],'png');
            saveas(FIGBERALL,[fname '.fig'],'fig');
            saveas(FIGBERALL,[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGBERALL);
        end
    end
    
    % ********************** ************************ *********************
    if (exist('hWB','var')&& ishandle(hWB))
        delete(hWB);
    end
catch ex
    if (exist('hWB','var')&& ishandle(hWB))
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

fprintf('--scrCSKPLLinear Done--\n');

















































