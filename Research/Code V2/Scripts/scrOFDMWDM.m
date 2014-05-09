% scrOFDMWDM
if exist('h','var') && ishandle(h)
    delete(h);
end
close all;
clearvars;
clc;

% FLAGS
fSTATION = 4;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = true;
fARCHIVE = false;
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = '1_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
NTX = 3; NRX = 3;
RMN = 627; RSD = 10; RSC = 1;
GMN = 530; GSD = 10; GSC = 1;
BMN = 470; BSD = 10; BSC = 1;
RES = 0.1;
sPSDTYP = 'Gaussian';
TOTALBITS = 5e4;
BERTH = 1e-3;   BERTHPL = 0.9*BERTH;
MODNSC = 64;

WBX = 50; WBY = 500; WBW = 275; WBH = 75;
WBTITLE = 'Running WDM OFDM Simulation...';

%% ranges
RNGCCT = 5500;                   LENCCT = numel(RNGCCT);
RNGCCTPL = RNGCCT;                      LENCCTPL = numel(RNGCCTPL);
RNGOFDMTYPES = {'dcoofdm';'acoofdm'};   LENOFDMTYPES = numel(RNGOFDMTYPES);
RNGOFDMOFST = [3.2;0];                  LENOFDMOFST = size(RNGOFDMOFST,2);
RNGMOD = [8;8^2];                       LENMOD  = size(RNGMOD,2);
RNGSNRMIN = 140; RNGSNRMAX = 300;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN;
BERRATIOS = [5 10 50 100 500]; DELTASNR = [0.05 0.1 2 4 8];

%% config
lkIl = 400;
lkPl = 1;       % plane for illumination
lkTx = 3;       % plane for transmitters
rmL = 4; rmW = 4; rmH = 4; % room L,W,H
rxX = 2; rxY = 2; rxZ = lkPl; % receiver location
flT = 1;
SIGBW = 40e6;

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\12 WDMOFDM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOFDMWDM.m';
        ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\';
        sPSDDIR = [ctMatDir sPSDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/12 WDMOFDM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOFDMWDM.m';
        ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/LEDPSD/';
        sPSDDIR = [ctMatDir sPSDTYP '/' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '/'];
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
        sPSDDIR = [ctMatDir sPSDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
    case 4
        ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\Pankil\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
        sPSDDIR = [ctMatDir sPSDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrOFDMWDM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes STRPREFIX 'datOFDMWDM' CHARIDXARCHIVE '.mat'];      % file to store workspace
RGBledmat = [sPSDDIR sprintf('res_%0.5f',RES) '.mat'];
if ~exist(ctDirRes,'dir')
    mkdir(ctDirRes);
end
%% PSDs

Rpsd = getSOG(RMN,RSD,RSC,lambdas);
Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);

Gpsd = getSOG(GMN,GSD,GSC,lambdas);
Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);

Bpsd = getSOG(BMN,BSD,BSC,lambdas);
Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);

Ambpsd = 5.8*ones(size(lambdas));        % Ambient PSD
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);   % Ambient Channel

%% variables
% TODO: NEED MORE RESOLUTION ON XYZ
if exist(RGBledmat,'file')
    load(RGBledmat,'RGB');
else
    RGB = cLEDrgb(RES,Rch,Gch,Bch);
    RGB.initialize();
    if ~exist(sPSDDIR,'dir')
        mkdir(sPSDDIR);
    end
    save(RGBledmat,'RGB');
end

% initialize receiver filter models
% TODO: May want to use 'Lorentzian' model
Rf = getSOG(RMN,RSD,RSC,lambdas);
Rf = flT*Rf/max(Rf);
Gf = getSOG(GMN,GSD,GSC,lambdas);
Gf = flT*Gf/max(Gf);
Bf = getSOG(BMN,BSD,BSC,lambdas);
Bf = flT*Bf/max(Bf);

% initialize receivers
Rrx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Rrx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
Grx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Grx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
Brx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Brx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);

%% logic

% result variable buffers
% BPS = zeros(LENOFDMTYPES,1);
% ERR = zeros(NTX,LENSNR,LENCCT,LENOFDMTYPES);
% BER = zeros(NTX,LENSNR,LENCCT,LENOFDMTYPES);
% XAVG = zeros(1,LENCCT);
% XAVGCH =  zeros(NTX,LENCCT);
ERR = [];
RNGSNRDB = [];
try
h = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
                'Visible','Off',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
set(h,'Position',[WBX WBY WBW WBH],'Visible','On');
setappdata(h,'canceling',0);
TOTALLOOPS = LENCCT*LENOFDMTYPES*RNGSNRLOOP;
for iT = 1:LENCCT
    % initialize transmitter
    [x,y] = planckXY(RNGCCT(iT));
    [~,R,G,B] = RGB.getPSD(x,y);
    
    % set room with new transmitter
    clear room;
    room = cRoom(rmL,rmW,rmH);        % create room and set size (4m x 4m x 4m)
    room.luminaire = cLuminaire([R G B],rmL/2,rmW/2,lkTx);
    room.luminaire.orientation = cOrientation(pi,0,0);
    locCntr = cLocation(rmL/2,rmW/2,lkPl);
    room.setIlluminance(locCntr,lkIl);        % set illuminance at receiver location
    % S = R+G+B; S = room.luminaire.color
    % R,G,B are the luminaire channels
    % lmFluxClr and rdFluxClr are the individual channel fluxes
    XAVGCH(:,iT) = [R.rdFlux;G.rdFlux;B.rdFlux];
    XAVG(iT) = sum(XAVGCH(:,iT),1);
    % calculate H red
    tHfs = room.getFreeSpaceGain(Rrx.location,Rrx.orientation,Rrx.rxFOV);
    [HRr,IAmR] = Rrx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HRg,~] = Rrx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HRb,~] = Rrx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H green
    tHfs = room.getFreeSpaceGain(Grx.location,Grx.orientation,Grx.rxFOV);
    [HGr,IAmG] = Grx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HGg,~] = Grx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HGb,~] = Grx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H blue
    tHfs = room.getFreeSpaceGain(Brx.location,Brx.orientation,Brx.rxFOV);
    [HBr,IAmB] = Brx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HBg,~] = Brx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HBb,~] = Brx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % construct channel matrix (does not include transmitted power)
    H = [HRr HRg HRb;...
        HGr HGg HGb;...
        HBr HBg HBb];
    H = H.*[1 0 0;0 1 0;0 0 1];
    
    % Loop through different configurations
    for iOf = 1:LENOFDMTYPES
        ofdmType = lower(RNGOFDMTYPES{iOf});
        switch lower(ofdmType)
            case 'acoofdm'
                d = MODNSC/4;
                STROFDM = 'ACO';
                OffsetDcoStddev = 0;
                OffsetAcoStddev = RNGOFDMOFST(iOf);
            case 'dcoofdm'
                d = MODNSC/2-1;
                STROFDM = 'DCO';
                OffsetDcoStddev = RNGOFDMOFST(iOf);
                OffsetAcoStddev = 0;
        end
        M = RNGMOD(iOf);
        SYMS = getQAMsyms(M);
        BPS(iOf) = d*log2(M);
        X = zeros(NTX,MODNSC);
        Y = zeros(NRX,MODNSC);
        
        % Monte-Carlo O-OFDM runs
        LOOPDONE = false;
%         for iSNR = 1:LENSNR
        iSNR = 1;
        while ~LOOPDONE
            if iSNR > 1
%                 if isequal(sum(BER(:,iSNR-1,iT,iOf) < BERTHPL),NTX-1)
%                     dber = 0;
%                 end
                dber = min(BER(BER(:,iSNR-1,iT,iOf) > BERTHPL,iSNR-1,iT,iOf));
%                 RNGSNRDB(iSNR,iT,iOf) = RNGSNRDB(iSNR-1,iT,iOf) + getDeltaSNR(dber,DSNR(iSNR-1,iT,iOf),-1e-2,DSNRMAX);;
                RNGSNRDB(iSNR,iT,iOf) = RNGSNRDB(iSNR-1,iT,iOf) + getDeltaSNR(BERTHPL,dber,BERRATIOS,DELTASNR);
            else
                RNGSNRDB(iSNR,iT,iOf) = RNGSNRMIN;
            end
            ERR(1:NTX,iSNR,iT,iOf) = 0;
%             SNRDB = RNGSNRDB(iSNR);
            BITSTART = 1;
            IBITSTX = 0;
            while(BITSTART <= TOTALBITS-BPS(iOf)+1)
                BITSTOP = BITSTART + BPS(iOf) - 1;
                sym_bits = randi([0 1],[NTX,BPS(iOf)]);
                
                % Generate symbols X
                for iTx = 1:NTX
                    sym_dec = bin2decMat(reshape(sym_bits(iTx,:),d,log2(M)))+1;
                    X(iTx,:) = genOFDMsignal(... % Variable Arguments to the function
                        'data',sym_dec,...
                        'OFDMtype',ofdmType,...
                        'N',MODNSC,...
                        'Symbols',SYMS,...
                        'OffsetDcoStddev', OffsetDcoStddev,...
                        'OffsetAcoStddev', OffsetAcoStddev)';
                end
                
                % Scale X to meet average illumination constraints
                Xmn = mean(X,2);
                XSCL = XAVGCH(:,iT)./Xmn;
                for iTx = 1:NTX
                    Xs(iTx,:) = X(iTx,:)*XSCL(iTx);
                end
                
                % Calculate noise
%                 SNR = power(10,SNRDB/10);
                SNR = power(10,RNGSNRDB(iSNR,iT,iOf)/10);
                SNRrt = sqrt(SNR);
                W = (XAVG(iT)/SNRrt)*randn(NRX,MODNSC);
                
                % Compute output
                Y = H*Xs + W;
                % Y = H*Xs;
                
                % Estimate X
                Xhs = H\Y;
                
                % Unscale Xh
                Xh = Xhs./repmat(XSCL,1,MODNSC);
                
                % Recover bit stream(s)
                sym_bits_h = -1*ones(size(sym_bits));
                for iTx = 1:NTX
                    sym_dec_h = decodeOFDMsignal(Xh(iTx,:),...
                        'OFDMtype',ofdmType,...
                        'N',MODNSC,...
                        'Symbols',SYMS);
                    sym_bits_mat_h = dec2binMat(sym_dec_h-1,log2(M));
                    sym_bits_h(iTx,:) = sym_bits_mat_h(:);
                    errs = biterr2(sym_bits_h(iTx,:),sym_bits(iTx,:));
                    ERR(iTx,iSNR,iT,iOf) = ERR(iTx,iSNR,iT,iOf) + errs;
%                     ERR(end) = ERR(end) + errs;
                end % Transmitters
                
                IBITSTX = IBITSTX + BPS(iOf);
                BITSTART = BITSTOP + 1;
            end % while loop
            BER(:,iSNR,iT,iOf) = ERR(:,iSNR,iT,iOf)/IBITSTX;
%             BER(:,iSNR,iT,iOf) = ERR(end)/IBITSTX;
            
            if iSNR>1
                DBER(:,iSNR,iT,iOf) = BER(:,iSNR,iT,iOf) - BER(:,iSNR-1,iT,iOf);
                DSNR(iSNR,iT,iOf) = RNGSNRDB(iSNR,iT,iOf) - RNGSNRDB(iSNR-1,iT,iOf);
            else
                DBER(:,iSNR,iT,iOf) = BER(:,iSNR,iT,iOf);
                DSNR(iSNR,iT,iOf) = 0;
            end
%             sum(BER(:,iSNR,iT,iOf) < BERTHPL)
            LOOPCOUNT = (iT*iOf-1)*RNGSNRLOOP + RNGSNRDB(iSNR,iT,iOf) - RNGSNRMIN;
            PROGRESS = LOOPCOUNT/TOTALLOOPS;
            waitbar(PROGRESS,h,sprintf('Simulation: %0.2f%% done...',PROGRESS*100));
            if(getappdata(h,'canceling'))
                delete(h);
                error('Simulation aborted');
            end
        
            if (RNGSNRDB(iSNR,iT,iOf) > RNGSNRMAX) || isequal(sum(BER(:,iSNR,iT,iOf) < BERTHPL),NTX)
                LOOPDONE = true;
                LOOPCOUNT = (iT*iOf)*RNGSNRLOOP;
                PROGRESS = LOOPCOUNT/TOTALLOOPS;
                waitbar(PROGRESS,h,sprintf('Simulation: %0.2f%% done...',PROGRESS*100));
                if(getappdata(h,'canceling'))
                    delete(h);
                    error('Simulation aborted');
                end
            end
            iSNR = iSNR + 1;
        end % SNR
        
    end % OFDM
end % CCT
% ERR = reshape(ERR,[NTX,
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end
delete(h);
catch ex
    delete(h);
    rethrow(ex);
end
scrOFDMWDMPL(ctFileVars);





































