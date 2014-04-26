% scrOSM_ofdm_NI
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
% STATION
switch fSTATION
    case 1
        ctFileCodeSrcIMG = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrSISOFDM_IMG.m';
    case 2
        ctFileCodeSrcIMG = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrSISOFDM_IMG.m';
    case 3
        ctFileCodeSrcIMG = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrSISOFDM_IMG.m';
    otherwise
        error('Station not defined');
end
run(ctFileCodeSrcIMG);
close all;
clearvars -except 'fSTATION';
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
% daxesfontsize = get(0,'DefaultAxesFontSize');
% set(0,'DefaultAxesFontSize',16);
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

% FLAGS
% fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fARCHIVE = true;
fFILTBLUE = false;

rand('seed',0); % seed for Random Number Generator
randn('seed',0); % seed for Random Number Generator

MODSSK = 1;
MODnSM = 2;
MODgSM = 3;
MODeSM = 4;

CHAROVERWRITE = '~';
STRPREFIX = '3_SISO_H_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
STRPREFIXIMG = '3_SIS_IMG_';
CHARIDXARCHIVEIMG = '';

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\11. SISO OFDM H\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrSISO_OFDM_H.m';
        ctDirResIMG = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\11. SISOFDM all\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/11. SISO OFDM H/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrSISO_OFDM_H.m';
        ctDirResIMG = '/home/pbutala/My Documents/MatlabResults/11. SISOFDM all/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\11. SISO OFDM H\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrSISO_OFDM_H.m';
        ctDirResIMG = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\11. SISOFDM all\\';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrSISO_OFDM_H' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes STRPREFIX 'scrSISO_OFDM_H' CHARIDXARCHIVE '.mat'];      % file to store workspace

% VARIABLES
%%% Setup
lkIl = 400;
lkPl = 1;       % plane for illumination
txa = 20e-2;     % transmitter side length
txD = 0.5;      % transmitter pitch
txm = 1;        % transmitter m
opFOV = 60*pi/180;   % optics FOV
rxal = 1e-3;    % pixel side length
rxD = 1e-3;
Nrxx = 1;
Nrxy = 1;
Nr = Nrxx*Nrxy;
rxZAT = cOrientation(0,0,0); % receiver orientation
%%% OSM
rngSNRdb = 150:0.25:350;
% rngSNRdb = [200 400];
lenSNRdb = length(rngSNRdb);
rngOfdmType = {'DCOOFDM','ACOOFDM'};
lenOfdmType = length(rngOfdmType);
rngOfstSDDcoS = 3.5;
rngOfstSDAcoS = 0;
rngOfstSDDcoM = 3.5;
rngOfstSDAcoM = 0.5;
% rngOfstSD = 0:0.5:5;
% rngOfstSD = 1:1:2;
lenOfstSD = length(rngOfstSDDcoS);

% rngMaco = [16 64];
% rngMdco = power(2,log2(sqrt(rngMaco)));
rngMacoS = 16384;
rngMdcoS = 128;
rngMacoM = 64;
rngMdcoM = 8;

lenM = length(rngMacoS);
Nsc = 64; % number of subcarriers

Ntx = 1;     % 1 transmitters
dtk = log2(Ntx);

% CONSTANTS
TOTALBITS = 5e4;
BERTH = 1e-3;
BITSTREAM = randi([0 1],[1,TOTALBITS]);
IDXBRK = 0;
STRSS = 'SISO';
% constants
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX; % wavelengths
s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;

Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);        % White PSD
txSPD = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);   % luminaire SPD
Ambpsd = 5.8e-2*ones(size(lambdas));        % Ambient PSD
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);   % Ambient Channel

% etc vars
clLC = {'k','b','r','g','m'};
lenLC = length(clLC);
clLS = {'-.','-','--',':'};
lenLS = length(clLS);
clMK = {'h','o','s','d','v','^','<','>','p'};
lenMK = length(clMK);

% delete old files
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end

if fFILTBLUE        % if blue filtering selected
    Bf = zeros(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
end

% initialize room
room = cRoom(4,4,4);        % create room and set size (4m x 4m x 4m)
locCntr = cLocation(room.L/2,room.W/2,1);
k = log2(Ntx);
Ntxx = power(2,ceil(k/2));
Ntxy = power(2,floor(k/2));
[txX, txY, txZ] = getGrid(Ntxx,Ntxy,1,txD,txD,2); % 2x2 array of transmitters
% [txX, txY, txZ] = getGrid(room.L,room.W,1,txD,txD,2,'Fill'); % 2x2 array of transmitter
txX = txX + room.L/2;       % add length location offset
txY = txY + room.W/2;       % add width location offset
txZ = txZ + 3;              % add height location offset
txLoc = cLocation(txX,txY,txZ); % tx location object
txOri = cOrientation(pi,pi/4,0);   % tx orientation
txSz = cSize(txa,txa,1e-3);     % tx size
for iL = 1:1:Ntx            % need in case transmitters have different output flux
    Wch(iL) = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);   % luminaire SPD
    aLum(iL) = cLuminaire(Wch(iL),cLocation(txX(iL),txY(iL),txZ(iL)));  % create luminaire object
    aLum(iL).order = txm;                       % set lamb order
    aLum(iL).orientation = txOri;               % set orientation
    aLum(iL).dimension = txSz;                  % set size
end
room.luminaire = aLum;      % place luminaires in room
room.setIlluminance(locCntr,lkIl);        % set illuminance at receiver location

% receiver setup
[rxX, rxY, rxZ] = getGrid(Nrxx,Nrxy,1,rxD,rxD,2); % 2x2 array of transmitters
rxX = rxX + room.L/2;       % add length location offset
rxY = rxY + room.W/2;       % add width location offset
rxZ = rxZ + 1;              % add height location offset
rxLoc = cLocation(rxX,rxY,rxZ);         %
Wrx = cSinglePixelReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z);
if fFILTBLUE
    Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);  % load blue filter
    Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1)); % receiver responsivity
end
Wrx.optics.fov = opFOV;                         % set FOV
Wrx.orientation = rxZAT;                        % set orientation

% calculate p-channel matrix
dStr = sprintf('Calculating H for Nt=%d x Np=%d',Ntx,Nr);
disp(dStr);

tHfs = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);  % free space channel gains
fsH = reshape(room.getFreeSpaceGain(rxLoc,cOrientation(0,0,0),pi/2),Wrx.rxCount,sum(room.lmCount(:)));
Hfs = reshape(tHfs,Wrx.rxCount,sum(room.lmCount(:)));   % reshape H to size=[Npx,Ntx]

% calculate channel matrix
color = room.luminaire(1).color;                    % get color of luminaires
txOri = room.luminaire(1).orientation;              % get tx orientations
[Hp,Iambpx] = ...
    Wrx.getSignal(txSPD,Hfs,Ambch);  % get p-channel matrix
% Hp = squeeze(tHp(1,:,:))';                        % reshape p-channel matrix
% if size(Hp,1) == 1
%     Hp = Hp';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OVERWRITE CHANNEL MATRIX
Hp = 0.1526e-7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BPSS = zeros(lenM,lenOfdmType);
BPSM = zeros(lenM,lenOfdmType);
% LEDLEN = zeros(lenM,lenOfdmType);
SYMLEN = zeros(lenM,lenOfdmType);
% bit_err_sym = zeros(lenSNRdb,lenM,lenOfdmType,lenOfstSD);
% bit_err_led = zeros(lenSNRdb,lenM,lenOfdmType,lenOfstSD);
bit_err = zeros(lenSNRdb,lenM,lenOfdmType,lenOfstSD);
figBER = zeros(lenOfstSD,lenM);
figBERall = zeros(lenOfstSD,lenM);
% figBERSplit = zeros(lenOfstSD,lenM);

OffsetDcoStddev = 0;
OffsetAcoStddev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iM = 1:lenM
    for iOfdm = 1:lenOfdmType
        ofdmType = rngOfdmType{iOfdm};
        switch lower(ofdmType)
            case 'acoofdm'
                Ms = rngMacoS(iM);
                Mm = rngMacoM(iM);
                d = Nsc/4;      % number of data carriers per ACOOFDM symbol
                STROFDM = 'ACO';
            case {'dcoofdm', 'dmt'}
                Ms = rngMdcoS(iM);
                Mm = rngMdcoM(iM);
                d = Nsc/2-1;      % number of data carriers per DCOOFDM/DMT symbol
                STROFDM = 'DCO';
        end
        dtm = log2(Ms);
        Syms = getQAMsyms(Ms);
        SYMLEN(iM,iOfdm) = d*dtm;
%         LEDLEN(iM,iOfdm) = Nsc*dtk;
%         BPS(iM,iOfdm) = SYMLEN(iM,iOfdm) + LEDLEN(iM,iOfdm);
        BPSS(iM,iOfdm) = SYMLEN(iM,iOfdm);
        BPSM(iM,iOfdm) = d*log2(Mm) + Nsc*2;
        for iOfst = 1:lenOfstSD
            OffsetDcoStddev = rngOfstSDDcoS(iOfst);
            OffsetAcoStddev = rngOfstSDAcoS(iOfst);
            STROFST = sprintf('OFST:%0.2f SD',rngOfstSDDcoS(iOfst));
            dStr = sprintf('%s, %s, BPS:%d, Nsc:%d, M:%d',STROFST,STROFDM,BPSS(iM,iOfdm),Nsc,Ms);
            disp(dStr);
            for idb = 1:lenSNRdb
                vSNRdb = rngSNRdb(idb);
                vSNR = power(10,vSNRdb/10);
                vSNRrt = sqrt(vSNR);
                err_sym = 0;
                err = 0;
                iBits = 0;
                
                while(iBits < TOTALBITS-BPSS(iM,iOfdm)+1);
                    SYMSTART = iBits+1;
                    SYMSTOP = SYMSTART+SYMLEN(iM,iOfdm)-1;
                    sym_bits = reshape(BITSTREAM(SYMSTART:SYMSTOP),d,dtm);
                    sym_index = bin2decMat(sym_bits)+1;
                    
%                     LEDSTART = SYMSTOP + 1;
%                     LEDSTOP = LEDSTART + LEDLEN(iM,iOfdm) - 1;
%                     led_bits = reshape(BITSTREAM(LEDSTART:LEDSTOP),Nsc,dtk);
%                     led_index = bin2decMat(led_bits)+1;
                    
                    tSig = genOFDMsignal(... % Variable Arguments to the function
                        'data',sym_index,...
                        'OFDMtype',ofdmType,...
                        'N',Nsc,...
                        'Symbols',Syms,...
                        'OffsetDcoStddev', OffsetDcoStddev,...
                        'OffsetAcoStddev', OffsetAcoStddev);
                    tSigMn = mean(tSig);
                    W = (tSigMn/vSNRrt)*randn(Nsc,Nr);
%                     X = zeros(Ntx,Nsc);
%                     X(sub2ind([Ntx,Nsc],led_index,(1:Nsc)')) = tSig;
                    X = tSig;
                    Y = Hp*X + W;
                    %             Y = Hp*X;
                    tSig_h = Hp\Y;
                    % decode LED index for SM
%                     [~, led_index_h] = max(tSig_h.*tSig_h,[],1);
%                     tSigI = sub2ind([Ntx,Nsc],led_index_h,1:Nsc);
%                     sym_index_h = decodeOFDMsignal(tSig_h(tSigI),...
%                         'OFDMtype',ofdmType,...
%                         'N',Nsc,...
%                         'Symbols',Syms);
                    sym_index_h = decodeOFDMsignal(tSig_h,...
                        'OFDMtype',ofdmType,...
                        'N',Nsc,...
                        'Symbols',Syms);
                    
                    sym_bits_h = dec2binMat(sym_index_h-1,dtm);
%                     led_bits_h = dec2binMat(led_index_h-1,dtk);
                    
                    err_sym = err_sym + biterr2(sym_bits_h,sym_bits);
%                     err_led = err_led + biterr2(led_bits_h,led_bits);
                    
%                     iBit_sym = iBit_sym + SYMLEN(iM,iOfdm);
%                     iBit_led = iBit_led + LEDLEN(iM,iOfdm);
                    iBits = iBits + BPSS(iM,iOfdm);
                end
%                 bit_err_sym(idb,iM,iOfdm,iOfst) = err_sym/iBit_sym;
%                 bit_err_led(idb,iM,iOfdm,iOfst) = err_led/iBit_led;
%                 err = err + err_sym + err_led;
                err = err + err_sym;
                bit_err(idb,iM,iOfdm,iOfst) = err/iBits;
                
%                 if isequal(bit_err_sym(idb,iM,iOfdm,iOfst),0);
%                         bit_err_sym(idb,iM,iOfdm,iOfst) = bit_err_sym(idb-1,iM,iOfdm,iOfst)/100;
%                 end
%                     
%                 if isequal(bit_err_led(idb,iM,iOfdm,iOfst),0);
%                         bit_err_led(idb,iM,iOfdm,iOfst) = bit_err_led(idb-1,iM,iOfdm,iOfst)/100;
%                 end
                    
                if bit_err(idb,iM,iOfdm,iOfst) < (BERTH/10)
                    if isequal(bit_err(idb,iM,iOfdm,iOfst),0);
                        bit_err(idb,iM,iOfdm,iOfst) = bit_err(idb-1,iM,iOfdm,iOfst)/100;
                    end
                    IDXBRK = IDXBRK + 1;
                    if(IDXBRK > 1)
                        bit_err(idb+1:end,iM,iOfdm,iOfst) = NaN;
                        IDXBRK = 0;
                        break;
                    end
                else
                    IDXBRK = 0;
                end
                dStr = sprintf('iter SNR:%0.2f BER:%0.5f',vSNRdb,bit_err(idb,iM,iOfdm,iOfst));
                disp(dStr);
            end % SNR
        end % Ofst
    end % Ofdm
end % M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOTDMIN = 4;
SNRD = zeros(lenOfstSD,lenM);
SNRA = zeros(lenOfstSD,lenM);
for iOfst = 1:lenOfstSD
    STROFST = sprintf('OfstSD_%0.2f',rngOfstSDDcoS(iOfst));
    for iM = 1:lenM
        STRMS = sprintf('_adM_%d_%d',rngMacoS(iM),rngMdcoS(iM));
        STRMM = sprintf('_adM_%d_%d',rngMacoM(iM),rngMdcoM(iM));
        clLgdall = {};
        figBER(iOfst,iM) = figure;
        set(gca,'YScale','log');
        hold all;
%         figBERSplit(iOfst,iM) = figure;
%         set(gca,'YScale','log');
%         hold all;
        fnameIMG = [ctDirResIMG STRPREFIXIMG STROFST STRMM '_ALL BER vs SNR' CHARIDXARCHIVEIMG];
        uiopen([fnameIMG '.fig'],1);
        figBERall(iOfst,iM) = gcf;
        [~,~,~,text_strings] = legend;
        for l=1:numel(text_strings)
            clLgdall{end+1} = text_strings{l};
        end
        
        clLgd = {};
%         clLgdSplit = {};
        iLC = 0;
        iLS = 0;
        iMK = 0;
        for iOfdm = 1:lenOfdmType
            ofdmType = rngOfdmType{iOfdm};
            iSnrTh = find(bit_err(:,iM,iOfdm,iOfst) <= BERTH,1,'first');
            switch lower(ofdmType)
                case 'acoofdm'
                    STROFDM = 'ACO';
                    STRMS = sprintf('M:%d',rngMacoS(iM));
                    if ~isempty(iSnrTh)
                        SNRA(iOfst,iM) = rngSNRdb(iSnrTh);
                    else
                        SNRA(iOfst,iM) = NaN;
                    end
                case {'dcoofdm', 'dmt'}
                    STROFDM = 'DCO';
                    STRMS = sprintf('M:%d',rngMdcoS(iM));
                    if ~isempty(iSnrTh)
                        SNRD(iOfst,iM) = rngSNRdb(iSnrTh);
                    else
                        SNRD(iOfst,iM) = NaN;
                    end
            end
            
            iLC = rem(iOfdm,lenLC)+1;
            iLS = rem(iM,lenLS)+1;
            iMK = rem(iMK+1,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            
            set(0,'CurrentFigure',figBER(iOfst,iM));
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            clLgd{end+1} = [STRSS ' ' STROFDM ' ' STRMS];
            clLgdall{end+1} = clLgd{end};
            legend(gca,clLgd,'Location','NorthEast');
%             xlabel('{P_{avg}^{tx}}/{N_{0}} (dB)');
            xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
            ylabel('BER');
            tStr = sprintf('BER vs SNR, N_{sc}:%d, Offset SISO DCO:%0.2f/ACO:%0.2f SD,\nDCO:%d/ACO:%d bits/sym',Nsc,rngOfstSDDcoS(iOfst),rngOfstSDAcoS(iOfst),BPSS(iM,1),BPSS(iM,2));
            title(tStr);
            grid on;
            axis([rngSNRdb(1)-150 rngSNRdb(end)-150 BERTH/5 1]);
            
            set(0,'CurrentFigure',figBERall(iOfst,iM));
            iLS = rem(3,lenLS)+1;
            iMK = rem(iMK+3,lenMK)+1;
            plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
            [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err(:,iM,iOfdm,iOfst)),PLOTDMIN);
%             semilogy(gca,rngSNRdb,bit_err(:,iM,iOfdm,iOfst),plStyle);
            semilogy(gca,Xc-150,power(10,Yc),plStyle);
            legend(gca,clLgdall,'Location','NorthEast');
            tStr = sprintf('BER vs SNR, N_{sc}:%d, Offset SISO DCO:%0.2f/ACO:%0.2f SD, SIS DCO:%0.2f/ACO:%0.2f SD\nSISO DCO:%d/ACO:%d bits/sym, SIS DCO:%d/ACO:%d bits/sym',...
                            Nsc,rngOfstSDDcoS(iOfst),rngOfstSDAcoS(iOfst),rngOfstSDDcoM(iOfst),rngOfstSDAcoM(iOfst),BPSS(iM,1),BPSS(iM,2),BPSM(iM,1),BPSM(iM,2));
            title(tStr);
%             set(0,'CurrentFigure',figBERSplit(iOfst,iM));
%             iLS = rem(2,lenLS)+1;
%             iMK = rem(iMK+1,lenMK)+1;
%             plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%             [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_sym(:,iM,iOfdm,iOfst)),PLOTDMIN);
% %             semilogy(gca,rngSNRdb,bit_err_sym(:,iM,iOfdm,iOfst),plStyle);
%             semilogy(gca,Xc-150,power(10,Yc),plStyle);
%             clLgdSplit{end+1} = [STRNI ' ' STROFDM ' ' STRM '- SYM'];
%             iLS = rem(3,lenLS)+1;
%             iMK = rem(iMK+1,lenMK)+1;
%             plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%             [Xc,Yc] = getCleanPoints(rngSNRdb,log10(bit_err_led(:,iM,iOfdm,iOfst)),PLOTDMIN);
% %             semilogy(gca,rngSNRdb,bit_err_led(:,iM,iOfdm,iOfst),plStyle);
%             semilogy(gca,Xc-150,power(10,Yc),plStyle);
%             clLgdSplit{end+1} = [STRNI ' ' STROFDM ' ' STRM '- TX'];
%             
%             legend(gca,clLgdSplit,'Location','NorthEast');
% %             xlabel('{P_{avg}^{tx}}/{N_{0}} (dB)');
%             xlabel('{SNR_{avg}^{tx}} - 150 (dB)');
%             ylabel('BER (Individual)');
%             tStr = sprintf('BER (Split) vs SNR, Offset: %0.2f SD\nN_{sc}:%d, DCO:%d/ACO:%d bits/sym',rngOfstSD(iOfst),Nsc,BPS(iM,1),BPS(iM,2));
%             title(tStr);
%             grid on;
%             axis([rngSNRdb(1)-150 rngSNRdb(end)-150 BERTH/5 1]);
        end % Ofdm
    end % M
end % Ofst

if lenOfstSD > 1
    figOFST = figure;
    clLgdOfst = {};
    hold all;
    iMK = 0;
    for iM = 1:lenM
        iLC = rem(1,lenLC)+1;
        iLS = rem(1,lenLS)+1;
        iMK = rem(iMK+1,lenMK)+1;
        plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
        plot(rngOfstSDDcoS,SNRD(:,iM),plStyle);
        STRMS = sprintf('M:%d',rngMdcoS(iM));
        clLgdOfst{end+1} = [STRSS ' ' 'DCO ' STRMS];
        
        iLC = rem(2,lenLC)+1;
        iMK = rem(iMK+1,lenMK)+1;
        plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
        plot(rngOfstSDDcoS,SNRA(:,iM),plStyle);
        STRMS = sprintf('M:%d',rngMacoS(iM));
        clLgdOfst{end+1} = [STRSS ' ' 'ACO ' STRMS];
    end
    xlabel('Offset (Std Dev)');
    yStr = sprintf('{SNR_{avg}^{tx}} - 150 (dB) (BER=10^{%d})',log10(BERTH));
    ylabel(yStr);
    tStr = sprintf('SNR vs Offset, Target BER = 10^{%d}',log10(BERTH));
    title(tStr);
    grid on;
    axis([rngOfstSDDcoS(1)-150 rngOfstSDDcoS(end)-150 rngSNRdb(1) rngSNRdb(end)]);
    legend(gca,clLgdOfst,'Location','NorthEast');
end

% draw setup
figSetup = figure;
room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
rotate3d on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
    fname = [ctDirRes STRPREFIX 'Setup' CHARIDXARCHIVE];
    saveas(figSetup,[fname '.png'],'png');
    saveas(figSetup,[fname '.fig'],'fig');
        saveas(figSetup,[fname '.eps'],'epsc');
    for iOfst = 1:lenOfstSD
        STROFST = sprintf('OfstSD_%0.2f',rngOfstSDDcoS(iOfst));
        for iM = 1:lenM
            STRMS = sprintf('_adM_%d_%d_%d_%d',rngMacoS(iM),rngMdcoS(iM),rngMacoM(iM),rngMdcoM(iM));
            f = figure(figBER(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRMS '_BER vs SNR' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
            f = figure(figBERall(iOfst,iM));
            fname = [ctDirRes STRPREFIX STROFST STRMS '_ALL BER vs SNR' CHARIDXARCHIVE];
            saveas(f,[fname '.png'],'png');
            saveas(f,[fname '.fig'],'fig');
            saveas(f,[fname '.eps'],'epsc');
%             f = figure(figBERSplit(iOfst,iM));
%             fname = [ctDirRes STRPREFIX STROFST STRM '_BER vs SNR Split' CHARIDXARCHIVE];
%             saveas(f,[fname '.png'],'png');
%             saveas(f,[fname '.fig'],'fig');
%             saveas(f,[fname '.eps'],'epsc');
        end
    end
    if lenOfstSD > 1
        fname = [ctDirRes STRPREFIX 'SNR vs Offset' CHARIDXARCHIVE];
        saveas(figOFST,[fname '.png'],'png');
        saveas(figOFST,[fname '.fig'],'fig');
        saveas(figOFST,[fname '.eps'],'epsc');
    end
end


% restore defaults
set(0,'DefaultLineMarkerSize',dlinems);
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
% set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);

disp(' ');
disp('--DONE--');






























