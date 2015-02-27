%scrMMPL
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
    RNGSNRMINPL = 0; RNGSNRMAXPL = 50;
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
    MKTYP = {'o','+','^','s','*','x','d','v','o','+','^','s','*','x','d','v'}; 
    MKCLR = {'g','y','b','r','g','b','y','r','r','b','y','g','g','y','b','r'};
    
    % Figure BER vs SNR config
    STRSNR = 'SNR';
    STRLSNR = 'LSNR';
    PLMMLCS = {'r';'g';'b';'r';'g';'b';'r';'g';'b'}; 
    PLMMLSS = {'--';'-.';':';'-.';':';'--';':';'--';'-.'}; 
    PLMMMKS = {'x';'+';'*';'s';'d';'h';'<';'v';'>'};
    PLMMLEN = numel(PLMMLCS);
    
    LOOPCOUNT = 0;
    TOTALLOOPS = LENMMCBC + (LENMMCBC>1)*LENMMCBC;
%     TOTALLOOPS = LENMMCBC*5 + (fSAVECHST==true)*sum(IDXCHST(:)) + (LENMMCBC>1)*LENMMCBC*2;
    TSTART = tic;
    
    % ********************** PLOT AND SAVE BER vs SNR *********************
    if LENMMCBC > 1
        sLGD = {}; hLGD = [];
        FIGBERALL = figure('Name',sprintf('%d-MM BER vs SNR',M),'NumberTitle',FIGTITLE);
        set(gca,'YScale','Log');
        hold on;
    end
    
    for iMMCBC=1:LENMMCBC
        PL1 = floor((iMMCBC-1)/PLMMLEN)+1;
        PL2 = rem((iMMCBC-1),PLMMLEN)+1;
        LSTL = [PLMMLCS{PL2} PLMMLSS{PL1}];
        MSTL = [PLMMLCS{PL2} PLMMMKS{PL1}];
        PSTL = [PLMMLCS{PL2} PLMMLSS{PL1} PLMMMKS{PL1}];
        FIGBER(iMMCBC) = figure('Name',sprintf('%d-MM CBs-%s BER vs SNR',M,mpC2Vstr{iMMCBC}),'NumberTitle',FIGTITLE);
        [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iMMCBC),log10(BER(:,iMMCBC)),PLOTDMIN);  % Get points well spaced out
        semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker
        hold on;
        semilogy(RNGSNROFST(:,iMMCBC),BER(:,iMMCBC),LSTL);          % Semilog AVG BER vs SNR
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        grid on;
        xlabel([STRSNR '(dB)']);
        ylabel('BER');
        title(sprintf('%d-MM BER vs SNR',M));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('%d-MM_CBs_%s_BERvsSNR',M,mpC2Vstr{iMMCBC}) CHARIDXARCHIVE];
            saveas(FIGBER(iMMCBC),[fname '.png'],'png');
            saveas(FIGBER(iMMCBC),[fname '.fig'],'fig');
%             saveas(FIGBER(iMMCBC),[fname '.eps'],'eps');
        end
        if fCLOSEALL
            close(FIGBER(iMMCBC));
        end
        % Update Wait bar
        LOOPCOUNT = LOOPCOUNT+1;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
        end
        
        if LENMMCBC > 1
            % BER vs SNR
            figure(FIGBERALL);
            semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker\\
            semilogy(RNGSNROFST(:,iMMCBC),BER(:,iMMCBC),LSTL);          % Semilog AVG BER vs SNR
            hLGD(end+1) = semilogy(nan,nan,PSTL);
            sLGD{end+1} = sprintf('CBs%s',mpC2Vstr{iMMCBC});
            
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
    if LENMMCBC > 1
        % snr
        figure(FIGBERALL);
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        legend(gca,hLGD,sLGD,'Location','NorthEast');
        grid on;
        xlabel([STRSNR '(dB)']);
        ylabel('BER');
        title(sprintf('%d-MM BER vs SNR',M));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('%d-MM_BERvsSNR',M) CHARIDXARCHIVE];
            saveas(FIGBERALL,[fname '.png'],'png');
            saveas(FIGBERALL,[fname '.fig'],'fig');
%             saveas(FIGBERALL,[fname '.eps'],'eps');
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
%     setpref('Internet','E_mail','pbutala@bu.edu');
%     setpref('Internet','SMTP_Server','smtp.bu.edu');
%     STREMAIL = ['Simulation ' STRPREFIX ' done with errors.'];
%     sendmail('pankil.butala@gmail.com',STREMAIL);
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
% setpref('Internet','E_mail','pbutala@bu.edu');
% setpref('Internet','SMTP_Server','smtp.bu.edu');
% STREMAIL = ['Simulation ' STRPREFIX ' done.'];
% sendmail('pankil.butala@gmail.com',STREMAIL);
% end
fprintf('--scrMMPL Done--\n');