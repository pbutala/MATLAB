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
    % Figure BER vs SNR config
    FIGBERYMIN = 0.9*BERTH; FIGBERYMAX = 1;
    
    LOOPCOUNT = 0;
    TOTALLOOPS = 2;
    
    % ********************** PLOT AND SAVE COLOR GAMUT*********************
    FIGGMT = figure('Name','RGB LED xy gamut','NumberTitle',FIGTITLE);
    RGB.obs.showGamut();
    hold on;
    scatter(RGB.xyz(1,:),RGB.xyz(2,:),'.','b');
    title(sprintf('RGB LED xy gamut\nR:%dnm G:%dnm B:%dnm Res:%0.2f',RMN,GMN,BMN,RES));
    if fSAVEALL
        fname = [ctDirRes STRPREFIX 'Gamut' CHARIDXARCHIVE];
%         saveas(FIGGMT,[fname '.png'],'png');
        saveas(FIGGMT,[fname '.fig'],'fig');
%         saveas(FIGGMT,[fname '.eps'],'epsc');
    end
    if fCLOSEALL
        close(FIGGMT);
    end
    % Update Wait bar
    LOOPCOUNT = LOOPCOUNT+1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
    
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
        fname = [ctDirRes STRPREFIX 'BERvsSNR' CHARIDXARCHIVE];
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
    waitbar(PROGRESS,hWB,sprintf('Results: %0.2f%% done...',PROGRESS*100));
    
    % ********************** ************************ *********************
    delete(hWB);
catch ex
    delete(hWB);
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