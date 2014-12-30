% scrCSKPL
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
    
    % LSNR COMPUTATION
    LDB = 20*log10(PTXAVG);                     % Get avg elec power in db
    LDBIDX = find(LDB==min(LDB),1,'first');  % Chose reference CBC (max elec power to generate K lumens)
    LDBOFST = LDB - LDB(LDBIDX);             % Compute Ldb = 20*log10(Pref/Pcbc)
    RNGLSNROFST = RNGSNROFST - repmat(LDBOFST',size(RNGSNROFST,1),1);   % compute LSNRs
    
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
    TOTALLOOPS = LENCBC*3 + (fSAVECHST==true)*sum(IDXCHST(:)) + (LENCBC>1)*LENCBC*2;
    TSTART = tic;
    % ********************** PLOT AND SAVE SPDs responsivities *********************
    for iCBC=1:LENCBC
        fCBC = RNGCBC(iCBC);
        FIGSPD = figure('Name',sprintf('Spectral Spreads CBC%d',fCBC),'NumberTitle',FIGTITLE);
        hold on;
        sLGD = {}; hLGD = [];
        for i=1:RGBLED(iCBC).NCLR
            p = RGBLED(iCBC).PSDs{i};
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
            fname = [ctDirData STRPREFIX sprintf('CBC%d_SpectralSpreads',fCBC) CHARIDXARCHIVE];
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
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
        end
    end
    
    % ********************** PLOT AND SAVE COLOR GAMUT*********************
    for iCBC=1:LENCBC
        fCBC = RNGCBC(iCBC);
        FIGGMT = figure('Name',sprintf('CBC%d RGB LED xy gamut',fCBC),'NumberTitle',FIGTITLE);
        RGBLED(iCBC).obs.showGamut();
        hold on;
        scatter(RGBLED(iCBC).xyz(1,:),RGBLED(iCBC).xyz(2,:),16,'.','m');
        [x,y,~] = RGBLED(iCBC).obs.getCoordinates(RGBLED(iCBC).PSDs{1});
        scatter(x,y,80,'x','r','LineWidth',2);
        [x,y,~] = RGBLED(iCBC).obs.getCoordinates(RGBLED(iCBC).PSDs{2});
        scatter(x,y,80,'x','g','LineWidth',2);
        [x,y,~] = RGBLED(iCBC).obs.getCoordinates(RGBLED(iCBC).PSDs{3});
        scatter(x,y,80,'x','b','LineWidth',2);
        
        title(sprintf('CBC%d RGB LED xy gamut\nR:%dnm G:%dnm B:%dnm Res:%0.2f',fCBC,RMN,GMN,BMN,RES));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('CBC%d_Gamut',fCBC) CHARIDXARCHIVE];
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
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
        end
    end
    % *********************** PLOT AND SAVE SYMBOLS ***********************
    if fSAVECHST
        for iCBC=1:LENCBC
            fCBC = RNGCBC(iCBC);
            for i=1:IDXCHST(iCBC)
                FileChnlSt = [ctFileChnlStPRE sprintf('_CBC%d_%d',fCBC,i) CHARIDXARCHIVE '.mat'];
                load(FileChnlSt);
                FIGCHST = figure('Name',sprintf('CBC%d Constellation',fCBC),'NumberTitle',FIGTITLE);
                
                switch fDECODER
                    case {2}
                        RGBLED(iCBC).obs.showGamut();
                end
                hold on;
                for j=1:size(CHST(iCBC).SYMS,2)
                    I = find(CHST(iCBC).TxIdx == j);
                    if ~isempty(I)
                        switch fDECODER
                            case {1,3}
                                scatter3(CHST(iCBC).RxSymEst(1,I),CHST(iCBC).RxSymEst(2,I),CHST(iCBC).RxSymEst(3,I),MKTYP{j},'MarkerEdgeColor',MKCLR{j},'linewidth',1);
                            case 2
                                scatter(CHST(iCBC).RxSymEst(1,I),CHST(iCBC).RxSymEst(2,I),MKTYP{j},'MarkerEdgeColor',MKCLR{j},'linewidth',1);
                        end
                    end
                end
                
                switch fDECODER
                    case 1
                        scatter3(CHST(iCBC).SYMS(1,:),CHST(iCBC).SYMS(2,:),CHST(iCBC).SYMS(3,:),80,'k','x','linewidth',2);
                        view(3);
                        rotate3d on;
                    case 2
                        scatter(CHST(iCBC).SYMS(1,:),CHST(iCBC).SYMS(2,:),80,'k','x','linewidth',2);
                    case 3
                        scatter3(CHST(iCBC).SYMS(1,:),CHST(iCBC).SYMS(2,:),CHST(iCBC).SYMS(3,:),80,'k','x','linewidth',2);
                        axis([0 1 0 1 0 1]);
                        view(3);
                        rotate3d on;
                end
                
                grid on;
                
                title(sprintf('CBC%d Signal constellation for SNR_{avg} = %0.2f(dB), BER = %0.1e',fCBC,CHST(iCBC).SNRdB,CHST(iCBC).BER));
                if fSAVEALL
                    fname = [ctDirData STRPREFIX sprintf('CBC%d_Constellation%d',fCBC,i) CHARIDXARCHIVE];
                    saveas(FIGCHST,[fname '.png'],'png');
                    saveas(FIGCHST,[fname '.fig'],'fig');
                    saveas(FIGCHST,[fname '.eps'],'epsc');
                end
                if fCLOSEALL
                    close(FIGCHST);
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
        end
    end
    % ********************** PLOT AND SAVE BER vs SNR and LSNR *********************
    if LENCBC > 1
        sLGD = {}; hLGD = [];
        FIGBERALL = figure('Name',sprintf('BER vs SNR'),'NumberTitle',FIGTITLE);
        set(gca,'YScale','Log');
        hold on;
        
        sLGDl = {}; hLGDl = [];
        FIGLSNR = figure('Name',sprintf('BER vs LSNR'),'NumberTitle',FIGTITLE);
        set(gca,'YScale','Log');
        hold on;
    end
    
    for iCBC=1:LENCBC
        fCBC = RNGCBC(iCBC);
        LSTL = [PLCBCLCS{iCBC} PLCBCLSS{iCBC}];
        MSTL = [PLCBCLCS{iCBC} PLCBCMKS{iCBC}];
        PSTL = [PLCBCLCS{iCBC} PLCBCLSS{iCBC} PLCBCMKS{iCBC}];
        FIGBER(iCBC) = figure('Name',sprintf('CBC%d BER vs SNR',fCBC),'NumberTitle',FIGTITLE);
        [Xp,Yp] = getCleanPoints(RNGSNROFST(:,iCBC),log10(BER(:,iCBC)),PLOTDMIN);  % Get points well spaced out
        semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker
        hold on;
        semilogy(RNGSNROFST(:,iCBC),BER(:,iCBC),LSTL);          % Semilog AVG BER vs SNR
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        grid on;
        xlabel([STRSNR '(dB)']);
        ylabel('BER');
        title(sprintf('CBC%d BER vs SNR',fCBC));
        if fSAVEALL
            fname = [ctDirData STRPREFIX sprintf('CBC%d_BERvsSNR',fCBC) CHARIDXARCHIVE];
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
            % BER vs SNR
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
            waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
            
            % BER vs LSNR
            figure(FIGLSNR);
            [Xp,Yp] = getCleanPoints(RNGLSNROFST(:,iCBC),log10(BER(:,iCBC)),PLOTDMIN);  % Get points well spaced out
            semilogy(Xp,power(10,Yp),MSTL);                         % Semilog AVG BER vs SNR Marker\\
            semilogy(RNGLSNROFST(:,iCBC),BER(:,iCBC),LSTL);          % Semilog AVG BER vs SNR
            hLGDl(end+1) = semilogy(nan,nan,PSTL);
            sLGDl{end+1} = sprintf('CBC%d (%0.2f lm/W_{avg})',fCBC,PAVGPERLM(iCBC).^(-1));
            
            % Update Wait bar (lsnr)
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
        % snr
        figure(FIGBERALL);
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        legend(gca,hLGD,sLGD,'Location','SouthEast');
        grid on;
        xlabel([STRSNR '(dB)']);
        ylabel('BER');
        title(sprintf('BER vs SNR'));
        if fSAVEALL
            fname = [ctDirData STRPREFIX 'CBCALL_BERvsSNR' CHARIDXARCHIVE];
            saveas(FIGBERALL,[fname '.png'],'png');
            saveas(FIGBERALL,[fname '.fig'],'fig');
            saveas(FIGBERALL,[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGBERALL);
        end
        
        % lsnr
        figure(FIGLSNR);
        axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
        legend(gca,hLGDl,sLGDl,'Location','NorthEast');
        grid on;
        xlabel([STRLSNR '(dB)']);
        ylabel('BER');
        title(sprintf('BER vs LSNR (reference CBC%d)',RNGCBC(LDBIDX)));
        if fSAVEALL
            fname = [ctDirData STRPREFIX 'CBCALL_BERvsLSNR' CHARIDXARCHIVE];
            saveas(FIGLSNR,[fname '.png'],'png');
            saveas(FIGLSNR,[fname '.fig'],'fig');
            saveas(FIGLSNR,[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGLSNR);
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