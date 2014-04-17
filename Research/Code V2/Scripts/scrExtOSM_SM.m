% scrExtOSM
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

rand('seed',0); % seed for Random Number Generator
randn('seed',0); % seed for Random Number Generator

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = 'range';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\5. Ext OSM SM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrExtOSM_SM.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/5. Ext OSM SM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrExtOSM_SM.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\5. Ext OSM SM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrExtOSM_SM.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrExtOSM_nSM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'ExtOSM_nSM' CHARIDXARCHIVE '.mat'];      % file to store workspace

% VARIABLES
rngSNRdb = 1:1:25;
rngMI = 1;
lenMI = length(rngMI);
% rngM = power(2,1:2);
% lenM = length(rngM);
% rngNt = power(2,1:2);
% lenNt = length(rngNt);
typOSM = OSMtype.NSMPAM;
BPS = 4;
vNt = 1:6;
vM = power(2,1:3);
[clNt clM clbpS] = bpS2NtM(typOSM,BPS,vNt,vM,false);
oNt = clNt{1};
oM = clM{1};
obpS = clbpS{1};
nNtM = length(oNt);
% CONSTANTS
TOTALBITS = 1e5;
BITSTREAM = randi([0 1],[1,TOTALBITS]);

% etc vars
clLC = {'k','b','r','g','m'};
lenLC = length(clLC);
clLS = {'--','-',':','-.'};
lenLS = length(clLS);
clMK = {'h','o','x','+','s','d','v','^','<','>','p'};
lenMK = length(clMK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. SSK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figBER = figure;
set(gca,'YScale','log');
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',6);
hold all;
% clLgd = cell(lenM*lenNt*lenMI,1);
clLgd = cell(nNtM,1);
iLgd = 1;
for iMI = 1:lenMI
    MI = rngMI(iMI);
    for iNtM = 1:nNtM
        Nt = oNt(iNtM);
        M = oM(iNtM);
        %     for iNt = 1:lenNt
        %         Nt = rngNt(iNt);
        %         for iM = 1:lenM
        %             M = rngM(iM);
        k = log2(Nt);
        m = log2(M);
        
        Xs = getPAMsyms(M,1,1);   % get PAM symbols
        clsyms = cell(1,Nt);
        [clsyms{:}] = deal(Xs(2:end)');
        tsyms = blkdiag(clsyms{:});
        tSYMBS = [zeros(Nt,1) full(tsyms)];
        nBPSYM = floor(log2(size(tSYMBS,2)));
        nSYM = power(2,nBPSYM);
        optSYM = getCodeOpt(tSYMBS);
        SYMBS = optSYM(:,1:nSYM);
        
        EXv = mean(SYMBS,2);                  % compare for equivalent average powers
        EX = mean(EXv);
        
        % COMPUTE BER VS SNR
        H = eye(Nt);        % TODO: introuce H from imaging receiver setup
        bit_err = zeros(length(rngSNRdb),1);
        
        for idb = 1:length(rngSNRdb)
            vSNRdb = rngSNRdb(idb);
            vSNR = power(10,vSNRdb/10);
            vSNRrt = sqrt(vSNR);
            err = 0;
            iBits = 1;
            while(iBits <= TOTALBITS-nBPSYM+1);
                bitsBin = BITSTREAM(iBits:iBits+nBPSYM-1);
                bitsDec = bin2decMat(bitsBin);
                iBits = iBits + nBPSYM;
                
                X = SYMBS(:,bitsDec + 1);
                
                W = (EX/vSNRrt)*randn(Nt,1);
                
                Y = H*X + W;
                % Y = H*X;
                
                DYSYM = repmat(Y,1,nSYM) - SYMBS;
                DISTYSYM = sum((DYSYM.^2),1);
                bitsDecH = find(DISTYSYM == min(DISTYSYM),1);
                bitsBinH = dec2binMat(bitsDecH-1,nBPSYM);
                
                err = err + biterr2(bitsBin,bitsBinH);
            end
            bit_err(idb) = err/(iBits-1);
        end
        iLC = rem(iLgd,lenLC)+1;
        iLS = rem(iLgd,lenLS)+1;
        iMK = rem(iLgd,lenMK)+1;
        plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
        semilogy(rngSNRdb,bit_err,plStyle);                % plot bit error vs snr
        %             clLgd{iLgd} = sprintf('N_{t}=%0.0f,M=%0.0f,MI=%0.2f',Nt,M,MI);
        clLgd{iLgd} = sprintf('nSM: N_{t}=%0.0f,M=%0.0f\n%0.0f b/sym',Nt,M,obpS(iNtM));
        iLgd = iLgd + 1;
    %     end % M
    % end % iNt
    end
end

hold off;
xlabel('SNR (dB)');
ylabel('BER');
grid on;
axis tight;
legend(clLgd,'Location','SouthWest');
set(0,'DefaultLineMarkerSize',dlinems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. GENERALIZED SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. EXTENDED SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delete old files
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end

if fSAVEALL
    fname = [ctDirRes ' nOSM BER vs SNR' CHARIDXARCHIVE];
    %     fname = [ctDirRes ' OSM BER vs SNR']; % Use this to archive
    saveas(gcf,[fname '.png'],'png');
    saveas(gcf,[fname '.fig'],'fig');
end

if fCLOSEFIGS
    close;
end

% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);







































