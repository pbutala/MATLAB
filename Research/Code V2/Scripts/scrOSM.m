% scrOSM
close all;
clearvars;
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');

% FLAGS
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fARCHIVE = true;

CHAROVERWRITE = '~';
if(fARCHIVE) 
    CHARIDXARCHIVE = 'mp1';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\5. OSM img\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOSM.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/5. OSM img/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOSM.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\5. OSM img\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOSM.m';
    otherwise
        error('Station not defined');
end

ctFileCodeDest = [ctDirRes 'scrOSM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'OSM' CHARIDXARCHIVE '.mat'];      % file to store workspace

% VARIABLES
Nt = 1:8;
rngM = power(2,1);

% COMPUTE #BITS/SYMBOL
lenM = length(rngM);
clLgd = cell(lenM*4+1,1);
figSE = figure;
hold all;
% subplot(2,1,2);
% AXES(2) = gca;
% hold all;

dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',6);

% etc vars
clLC = {'k','b','r','g','m'};
lenLC = length(clLC);
clLS = {'--','-',':','-.'};
lenLS = length(clLS);
clMK = {'h','o','x','+','s','d','v','^','<','>','p'};
lenMK = length(clMK);
iLgd = 1;

for idM = 1:lenM
    M = rngM(idM);
    % SE for native SM (Ideal and 0:M-1 PAM)
    senSM = floor(log2(Nt*M));
    senSMp = floor(log2(1+Nt*(M-1)));
    % SE for generalized SM (Ideal and 0:M-1 PAM)
    segSM = floor(Nt + log2(M));
    segSMp = floor(log2(1+(power(2,Nt)-1)*(M-1)));
    % SE for extended SM (Ideal and 0:M-1 PAM)
    seeSYM = ones(length(Nt),1);
    seeSYMp = ones(length(Nt),1);
    for iNt = 1:length(Nt)
        nt = Nt(iNt);
        for k=1:nt
            seeSYM(iNt) = seeSYM(iNt) + nchoosek(nt,k)*power(M,k);
            seeSYMp(iNt) = seeSYMp(iNt) + nchoosek(nt,k)*power(M-1,k);
        end
    end
    seeSM = floor(log2(seeSYM));
    seeSMp = floor(log2(seeSYMp));
    
    iMK = rem(iLgd,lenMK)+1;
    
    stairs(Nt,seeSM,['g-' clMK{iMK}]);
    clLgd{iLgd} = sprintf('e-SM %0.0f-ary',M);
    stairs(Nt,seeSMp,['g:' clMK{iMK+1}]);
    clLgd{iLgd+1} = sprintf('e-SM %0.0f-PAM',M);
    stairs(Nt,segSM,['r-' clMK{iMK+2}]);
    clLgd{iLgd+2} = sprintf('g-SM %0.0f-ary',M);
    stairs(Nt,segSMp,['r:' clMK{iMK+3}]);
    clLgd{iLgd+3} = sprintf('g-SM %0.0f-PAM',M);
    stairs(Nt,senSM,['b-' clMK{iMK+4}]);
    clLgd{iLgd+4} = sprintf('n-SM %0.0f-ary',M);
    stairs(Nt,senSMp,['b:' clMK{iMK+5}]);
    clLgd{iLgd+5} = sprintf('n-SM %0.0f-PAM',M);
    
    iLgd = iLgd+6;
    iMK = iMK+6;
end

% SE for SSK
seSM = floor(log2(Nt));
stairs(Nt,seSM,'k-.h');
clLgd{iLgd} = 'SSK';

legend(clLgd{:},'Location','NorthWest');
xlabel('Number of transmitters (N_{t})');
ylabel('#-bits / symbol');
set(0,'DefaultLineMarkerSize',dlinems);
grid on;
axis([min(Nt),max(Nt),0,25]);

if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end

if fSAVEALL
    figure(figSE);
    fname = [ctDirRes ' OSMspeff' CHARIDXARCHIVE];
%     fname = [ctDirRes ' OSMspeff']; % Use this to archive
    saveas(gcf,[fname '.png'],'png');
    saveas(gcf,[fname '.fig'],'fig');
end
if fCLOSEFIGS
    close;
end

bpS2NtM(1:4,1,Nt,rngM,false);

% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);







































