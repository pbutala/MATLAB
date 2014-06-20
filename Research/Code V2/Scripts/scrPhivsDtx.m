% scrPhivsDtx
close all;
clearvars;
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',10);
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
FIGTITLE = 'Off';

try
    % FLAGS
    fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
    fSAVEALL = true;
    fCLOSEALL = true;
    fARCHIVE = true;
    
    CHAROVERWRITE = '~';
    STRPREFIX = '1_';
    if(fARCHIVE)
        CHARIDXARCHIVE = '';           % ARCHIVE INDEX
    else
        CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
    end
    
    % STATION
    switch fSTATION
        % Results directory; Spource file; LED table dir;
        case 1
            ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\13. PHIvsDtx\';
            ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrPhivsDtx.m';
        case 2
            ctDirRes = '/home/pbutala/My Documents/MatlabResults/13. PHIvsDtx/';
            ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrPhivsDtx.m';
        case 3
            ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\13. PHIvsDtx\\';
            ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrPhivsDtx.m';
        case 4
            ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\13. PHIvsDtx\\';
            ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrPhivsDtx.m';
        otherwise
            error('Station not defined');
    end
    ctFileCodeDest = [ctDirRes STRPREFIX 'scrPhivsDtx' CHARIDXARCHIVE '.m']; % Script copy name
    ctFileVars = [ctDirRes STRPREFIX 'datPhivsDtx' CHARIDXARCHIVE '.mat'];   % Data file name
    if ~exist(ctDirRes,'dir')   % if data directory does NOT exist
        mkdir(ctDirRes);        % create data dir
    end
    
    %% ranges
    RNGDtx = 0.05:0.05:5;                  LENDtx = numel(RNGDtx);             % Range of transmitter pitch
    RNGDNM = 0.5:0.5:4;                       LENDNM = numel(RNGDNM);             % Range of ditances between transmitter and receiver location planes
    
    RNGPHIs = zeros(LENDtx,LENDNM);
    %% logic
    for iDnm = 1:LENDNM
        dNM = RNGDNM(iDnm);
        for iDtx = 1:LENDtx
            Dtx = RNGDtx(iDtx);
            dPL = Dtx/sqrt(2);
            RNGPHIs(iDtx,iDnm) = atan(dPL/dNM);
        end
    end
    
    %% PLOT and SAVE
    PLC = {'b';'r';'g';'m'};       LENPLC = numel(PLC);
    PLS = {':';'-';'--';'-.'};         LENPLS = numel(PLS);
    LGD = {};
    
    if fSAVEALL                                                                 % SAVE script and data
        delete([ctDirRes '*' CHAROVERWRITE '.*']);
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
    
    FIGPHIALL = figure('Name',sprintf('Steer angle vs transmitter pitch'),'NumberTitle',FIGTITLE);
    axis([0 max(RNGDtx) 0 90]);
    xlabel('Transmitter Grid Pitch (m)');
    ylabel('Minimum Steer Half Angle (degrees)');
    title(sprintf('Steer angle vs transmitter pitch'));
    hold on;
    grid on;
    for iDnm = 1:LENDNM
        dNM = RNGDNM(iDnm);
        FIGPHI(iDnm) = figure('Name',sprintf('Steer angle vs transmitter pitch for receiver at %0.1fm distance',dNM),'NumberTitle',FIGTITLE);
        plot(RNGDtx,RNGPHIs(:,iDnm)*180/pi);
        axis([0 max(RNGDtx) 0 90]);
        grid on;
        xlabel('Transmitter Grid Pitch (m)');
        ylabel('Minimum Steer Half Angle (degrees)');
        title(sprintf('Steer angle vs transmitter pitch for receiver at %0.1fm distance',dNM));
        
        if fSAVEALL
            f = figure(FIGPHI(iDnm));
            fl = sprintf('PHIvsDTX_%0.1fmDNMrx',dNM);
            fname = [ctDirRes STRPREFIX fl CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
        end
        if fCLOSEALL
            close(FIGPHI(iDnm));
        end
       
        figure(FIGPHIALL);
        LC = PLC{rem(iDnm,LENPLC) + 1};
        LS = PLS{rem(floor((iDnm-1)/LENPLC)+1,LENPLS)+1};
        plot(RNGDtx,RNGPHIs(:,iDnm)*180/pi,[LC LS]);
        LGD{end+1} = sprintf('d_{rx}=%0.1fm',dNM); 
    end
    figure(FIGPHIALL);
    legend(LGD,'Location','SouthEast','FontSize',8);
    if fSAVEALL
        f = figure(FIGPHIALL);
        fl = sprintf('PHIvsDTX');
        fname = [ctDirRes STRPREFIX fl CHARIDXARCHIVE];
        saveas(f,[fname '.png'],'png');
        saveas(f,[fname '.fig'],'fig');
        saveas(f,[fname '.eps'],'epsc');
    end
    if fCLOSEALL
        close(FIGPHIALL);
    end
catch ex
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

% end