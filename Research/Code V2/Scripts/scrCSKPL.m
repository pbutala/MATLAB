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
    hWB = waitbar(0,'Plotting Results: 0.00% done','Name',WBTITLE,'Visible','Off');
    set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
    % PLOT Configs
    PLOTDMIN = 5;
    RNGSNROFST = RNGSNRDB - SNROFST;
    FIGBERXMIN = RNGSNRMIN-SNROFST; FIGBERXMAX = RNGSNRMAX-SNROFST;
    
    PLACOLC = 'm'; PLACOLS = '-'; PLACOMK = 'o';
    PLDCOLC = 'c'; PLDCOLS = '-'; PLDCOMK = 'd';
    PLTXLCS = {'r';'g';'b'}; PLTXLSS = {'--';'-.';':'}; PLTXMKS = {'>';'s';'*'};
    PLNSCMKS = {'x';'h';'^';'+';'v';'*';'<';'p'};
    MKTYP = {'o','+','^','s'}; 
    MKCLR = {'g','y','b','r'};
    % MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};
    % Figure BER vs SNR config
    FIGBERYMIN = 0.9*BERTH; FIGBERYMAX = 1;
    
    LOOPCOUNT = 0;
    TOTALLOOPS = 3 + IDXCHST;
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
    waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));
    
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
    waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));
    
    % *********************** PLOT AND SAVE SYMBOLS ***********************
    for i=1:IDXCHST
        FileChnlSt = [ctFileChnlStPRE sprintf('%d',i) CHARIDXARCHIVE '.mat'];
        load(FileChnlSt);
        FIGCHST = figure('Name','Constellation','NumberTitle',FIGTITLE);
        
        switch fDECODER
            case {2}
            RGB.obs.showGamut();
        end
        hold on;
        for j=1:size(CHST.SYMS,2)
            I = find(CHST.TxIdx == j);
            if ~isempty(I)
                switch fDECODER
                    case {1,3}
                        scatter3(CHST.RxSymEst(1,I),CHST.RxSymEst(2,I),CHST.RxSymEst(3,I),MKTYP{j},'MarkerEdgeColor',MKCLR{j},'linewidth',1);
                    case 2
                        scatter(CHST.RxSymEst(1,I),CHST.RxSymEst(2,I),MKTYP{j},'MarkerEdgeColor',MKCLR{j},'linewidth',1);
                end
            end
        end
        
        switch fDECODER
            case 1
                scatter3(CHST.SYMS(1,:),CHST.SYMS(2,:),CHST.SYMS(3,:),80,'k','x','linewidth',2);
                view(3);
                rotate3d on;
            case 2
                scatter(CHST.SYMS(1,:),CHST.SYMS(2,:),80,'k','x','linewidth',2);
            case 3
                scatter3(CHST.SYMS(1,:),CHST.SYMS(2,:),CHST.SYMS(3,:),80,'k','x','linewidth',2);
                axis([0 1 0 1 0 1]);
                view(3);
                rotate3d on;
        end
        
        grid on;
        
        title(sprintf('Signal constellation for SNR_{avg}^{tx} = %0.2f(dB), BER = %0.1e',CHST.SNRdB,CHST.BER));
        if fSAVEALL
            fname = [ctDirData STRPREFIX 'Constellation' sprintf('%d',i) CHARIDXARCHIVE];
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
        waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));
    end
    
    % ********************** PLOT AND SAVE BER vs SNR *********************
    FIGBER = figure;
    [Xp,Yp] = getCleanPoints(RNGSNROFST,log10(BER),PLOTDMIN);  % Get points well spaced out
    semilogy(Xp,power(10,Yp),'bo');                          % Semilog AVG BER vs SNR Marker
    hold on;
    semilogy(RNGSNROFST,BER,'b-');  % Semilog AVG BER vs SNR
    axis([FIGBERXMIN FIGBERXMAX FIGBERYMIN FIGBERYMAX]);
    grid on;
    xlabel(sprintf('SNR^{tx}_{avg}(dB)'));
    ylabel('BER');
    if fSAVEALL
        fname = [ctDirData STRPREFIX 'BERvsSNR' CHARIDXARCHIVE];
        saveas(FIGBER,[fname '.png'],'png');
        saveas(FIGBER,[fname '.fig'],'fig');
        saveas(FIGBER,[fname '.eps'],'epsc');
    end
    if fCLOSEALL
        close(FIGBER);
    end
    % Update Wait bar
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,hWB,sprintf('Plotting Results: %0.2f%% done...',PROGRESS*100));
    
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