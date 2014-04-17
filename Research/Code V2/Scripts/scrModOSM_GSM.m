% scrModOSM_GSM
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
fARCHIVE = false;

rand('seed',0); % seed for Random Number Generator
randn('seed',0); % seed for Random Number Generator

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = '2';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\6. Mod OSM GSM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrModOSM_GSM.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/6. Mod OSM GSM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrModOSM_GSM.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\6. Mod OSM GSM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrModOSM_GSM.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrModOSM_GSM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'ModOSM_GSM' CHARIDXARCHIVE '.mat'];      % file to store workspace

% VARIABLES
rngSNRdb = 1:1:25;
rngM = power(2,1);
lenM = length(rngM);
% rngMI = [1 0.95];
rngMI = 1;
lenMI = length(rngMI);
rngNt = power(2,1);
lenNt = length(rngNt);
% CONSTANTS
TOTALBITS = 1e1;
BITSTREAM = randi([0 1],[1,TOTALBITS]);

% seeSM = Nt + (Nt/2)*m;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. GENERALIZED SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figBER = figure;
set(gca,'YScale','log');
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',6);
hold all;
clLgd = cell(lenM*lenNt*lenMI,1);
iLgd = 1;
for iMI = 1:lenMI
    MI = rngMI(iMI);
    for iNt = 1:lenNt
        Nt = rngNt(iNt);
        for iM = 1:lenM
            M = rngM(iM);
%             k = log2(Nt);
            k = Nt; % for g-OSM
            m = log2(M);
            
            segSM = Nt + m;
            mk = segSM;
            
            % create all symbols
            Xs = getPAMsyms(M,2*MI/(M-1),1);
            EX2 = var(Xs);   % NOTE THIS IS DIFF THAN CONSIDERED IN RF
            EX = sqrt(EX2);     % just to calculate RELATIVE NOISE POWER
            % Xs = getPAMsyms(M,2*MI/(M-1),0);
            
            SYMBSNT = [zeros(Nt,1) reshape(repmat(dec2binMat(1:(power(2,Nt)-1)).',M-1,1),Nt,(M-1)*(power(2,Nt)-1))];
            SYMBSMOD = [zeros(Nt,1) repmat(repmat(Xs(2:end).',Nt,1),1,power(2,Nt)-1)];
            SYMBS = SYMBSNT.*SYMBSMOD;
            
            % COMPUTE BER VS SNR
            H = eye(Nt);
            bit_err = zeros(length(rngSNRdb),1);
            clear dataBitsBin dataBitsDec txBitsBin txBitsDec;
            clear dataBitsBinh dataBitsDech txBitsBinh txBitsDech;
            for idb = 1:length(rngSNRdb)
                vSNRdb = rngSNRdb(idb);
                vSNR = power(10,vSNRdb/10);
                vSNRrt = sqrt(vSNR);
                err = 0;
                iBits = 1;
                while(iBits <= TOTALBITS-mk+1);
                    dataBitsBin = BITSTREAM(iBits:iBits+m-1);
                    dataBitsDec = bin2decMat(dataBitsBin);
                    iBits = iBits + m;
                    if dataBitsDec == 0
                        X = SYMBS(:,1);
                        txBitsBin = [];
                    else
                        txBitsBin = BITSTREAM(iBits:iBits+k-1);
                        txBitsDec = bin2decMat(txBitsBin);
                        iBits = iBits + k;
                        
                        frmBitsBin = [dataBitsBin txBitsBin];
                        X = SYMBS(:,txBitsDec*M + dataBitsDec + 1);
                    end
                    W = (EX/vSNRrt)*randn(Nt,1);
                    
                    Y = H*X + W;
                    
                    DYSYM = repmat(Y,1,size(SYMBS,2)) - SYMBS;
                    DISTYSYM = sum((DYSYM.^2),1);
                    Xhidx = find(DISTYSYM == min(DISTYSYM),1);
                    
                    txBitsDech = floor((Xhidx-1)/M);
                    txBitsBinh = dec2binMat(txBitsDech,k);
                    
                    dataBitsDech =  rem(Xhidx-1,M);
                    dataBitsBinh = dec2binMat(dataBitsDech,m);
                    
                    estBits = [dataBitsBinh txBitsBinh];
                    
                    err = err + biterr2(frmBitsBin,estBits);
                end
                bit_err(idb) = err/(iBits-1);
            end
            iLC = rem(iLgd,lenLC)+1;
            iLS = rem(iLgd,lenLS)+1;
            iMK = rem(iLgd,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            semilogy(rngSNRdb,bit_err,plStyle);                % plot bit error vs snr
%             clLgd{iLgd} = sprintf('N_{t}=%0.0f,M=%0.0f,MI=%0.2f',Nt,M,MI);
            clLgd{iLgd} = sprintf('gSM: N_{t}=%0.0f,M=%0.0f',Nt,M);
            iLgd = iLgd + 1;
        end
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
    fname = [ctDirRes ' gOSM BER vs SNR' CHARIDXARCHIVE];
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







































