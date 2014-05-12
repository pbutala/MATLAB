% scrOFDMWDM
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% FLAGS
fSTATION = 4;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = true;
fARCHIVE = true;
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = '3_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
NTX = 3; NRX = 3;
RMN = 627; RSD = 10; RSC = 1;               % Mean, SD and scale to generate SPD of Red led
GMN = 530; GSD = 10; GSC = 1;               % Mean, SD and scale to generate SPD of Green led
BMN = 470; BSD = 10; BSC = 1;               % Mean, SD and scale to generate SPD of Blue led
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
sSPDTYP = 'Gaussian';                       % LED SPD model
WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running WDM OFDM Simulation...'; % Wait Box title

%% ranges
RNGCCT = [3000:500:6000];                         LENCCT = numel(RNGCCT);              % CCT 
RNGCCTPL = RNGCCT;                     LENCCTPL = numel(RNGCCTPL);          % CCTs to plot
RNGOFDMTYPES = {'dcoofdm';'acoofdm'};  LENOFDMTYPES = numel(RNGOFDMTYPES);  % OFDM types
RNGOFDMOFST = [3.2 0];                 LENOFDMOFST = size(RNGOFDMOFST,2);   % OFDM offsets
RNGMOD = [8 8^2];                      LENMOD  = size(RNGMOD,2);            % Subcarrier modulation order
RNGMODNSC = power(2,6);           LENMODNSC = numel(RNGMODNSC);        % Number of subcarriers
RNGSNRMIN = 140; RNGSNRMAX = 300;                                           % SNR ranges 
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000]; DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
% DELTASNR = [5 5 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR

TOTALBITS = 5e4;                            % Total bit for transmtter to simulate
BERTH = 1e-3;   BERTHMIN = 0.5*BERTH;       % BER thresholds; 

%% config
LKILL = 400;                    % Illumination Level
LKPL = 1;                       % plane for illumination
LKTX = 3;                       % plane for transmitters
RML = 4; RMW = 4; RMH = 4;      % room L,W,H
RXX = 2; RXY = 2; RXZ = LKPL;   % receiver location
FLT = 1;                        % Filter Transmission

% STATION
switch fSTATION
    % Results directory; Spource file; LED table dir; 
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\12 WDMOFDM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOFDMWDM.m';
        ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\';
        sPSDDIR = [ctMatDir sSPDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/12 WDMOFDM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOFDMWDM.m';
        ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/LEDPSD/';
        sPSDDIR = [ctMatDir sSPDTYP '/' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '/'];
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
        sPSDDIR = [ctMatDir sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
    case 4
        ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\Pankil\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
        sPSDDIR = [ctMatDir sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrOFDMWDM' CHARIDXARCHIVE '.m']; % Script copy name
ctFileVars = [ctDirRes STRPREFIX 'datOFDMWDM' CHARIDXARCHIVE '.mat'];   % Data file name
RGBledmat = [sPSDDIR sprintf('res_%0.5f',RES) '.mat'];                  % LED table mat-file  
if ~exist(ctDirRes,'dir')   % if data directory does NOT exist
    mkdir(ctDirRes);        % create data dir
end

%% PSDs
Rpsd = getSOG(RMN,RSD,RSC,lambdas);                 
Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);   % Normalized RED SPD

Gpsd = getSOG(GMN,GSD,GSC,lambdas);
Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);   % Normalized GREEN SPD

Bpsd = getSOG(BMN,BSD,BSC,lambdas);
Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);   % Normalized BLUE SPD

Ambpsd = 5.8*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);   % AMBIENT SPD

%% variables
% TODO: NEED MORE RESOLUTION ON XYZ
if exist(RGBledmat,'file')              % If LED characterization table exists
    load(RGBledmat,'RGB');              % Load RGB led
else
    RGB = cLEDrgb(RES,Rch,Gch,Bch);     % Create RGB led
    RGB.initialize();                   % Initialize led
    if ~exist(sPSDDIR,'dir')
        mkdir(sPSDDIR);                 % Create directory to store characterization data
    end
    save(RGBledmat,'RGB');              % Save LED characterization
end

% initialize receiver filter models
% TODO: May want to use 'Lorentzian' model
Rf = getSOG(RMN,RSD,RSC,lambdas);       % Red filter
Rf = FLT*Rf/max(Rf);
Gf = getSOG(GMN,GSD,GSC,lambdas);       % Green filter
Gf = FLT*Gf/max(Gf);
Bf = getSOG(BMN,BSD,BSC,lambdas);       % Blue filter
Bf = FLT*Bf/max(Bf);

% initialize receivers
Rrx = cSinglePixelReceiverWhiteReflection(RXX,RXY,RXZ);
Rrx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
Grx = cSinglePixelReceiverWhiteReflection(RXX,RXY,RXZ);
Grx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
Brx = cSinglePixelReceiverWhiteReflection(RXX,RXY,RXZ);
Brx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);

%% logic

% result variable buffers
ERR = []; RNGSNRDB = [];
try
% Wait Bar to show progress
hWB = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
    'Visible','Off',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
setappdata(hWB,'canceling',0);
LOOPCOUNT = 0;
TOTALLOOPS = LENCCT*LENMODNSC*LENOFDMTYPES*RNGSNRLOOP;
for iT = 1:LENCCT                                                           % LOOP START CCT
    % initialize transmitter
    [x,y] = planckXY(RNGCCT(iT));   % Get x,y from CCT
    [~,R,G,B] = RGB.getPSD(x,y);    % Get transmitter(s) SPD for set x,y
    % R,G,B are the luminaire channels
    
    % set room with new transmitter
    clear room;
    room = cRoom(RML,RMW,RMH);                                              % create room and set size 
    room.luminaire = cLuminaire([R G B],RML/2,RMW/2,LKTX);                  % Add luminaire to room
    room.luminaire.orientation = cOrientation(pi,0,0);                      % Set luminaire orientation
    
    locCntr = cLocation(RML/2,RMW/2,LKPL);                                  % Select center of room at plane (1m)
    room.setIlluminance(locCntr,LKILL);                                     % set illuminance
    
    % lmFluxClr and rdFluxClr are the individual channel fluxes
    XAVGCH(:,iT) = [R.rdFlux;G.rdFlux;B.rdFlux];                            % AVERAGE FLUX PER COLOR
    XAVG(iT) = sum(XAVGCH(:,iT),1);                                         % TOTAL AVERAGE FLUX
    
    % calculate H red
    tHfs = room.getFreeSpaceGain(Rrx.location,Rrx.orientation,Rrx.rxFOV);   % RED channel gains
    [HRr,~] = Rrx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HRg,~] = Rrx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HRb,~] = Rrx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H green
    tHfs = room.getFreeSpaceGain(Grx.location,Grx.orientation,Grx.rxFOV);   % GREEN channel gains
    [HGr,~] = Grx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HGg,~] = Grx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HGb,~] = Grx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H blue
    tHfs = room.getFreeSpaceGain(Brx.location,Brx.orientation,Brx.rxFOV);   % BLUE channel gains
    [HBr,~] = Brx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HBg,~] = Brx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HBb,~] = Brx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % construct channel matrix (does not include transmitted power)
    H = [HRr HRg HRb;...
        HGr HGg HGb;...
        HBr HBg HBb];                                                       % Channel matrix

    % Loop through different configurations
    for iNSC = 1:LENMODNSC                                                  % LOOP START MODNSC
        MODNSC = RNGMODNSC(iNSC);
        for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
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
            SYMS = getQAMsyms(M);                                           % Get QAM symbols
            BPS(iNSC,iOf) = d*log2(M);                                      % Calculate Bits Per Symbol
            
            X = zeros(NTX,MODNSC);
            Y = zeros(NRX,MODNSC);
            clear Xs Xhs Xh W;
            
            LOOPDONE = false; iSNR = 1;
            while ~LOOPDONE                                                 % LOOP START SNR (dynamically select next SNR)
                if iSNR > 1
                    dber = min(BER(BER(:,iSNR-1,iT,iNSC,iOf) > BERTHMIN,iSNR-1,iT,iNSC,iOf));   % Smallest BER from all those above the BERTH limit
                    % dynamically select next SNR based on how close the current BER is to threshold
                    RNGSNRDB(iSNR,iT,iNSC,iOf) = RNGSNRDB(iSNR-1,iT,iNSC,iOf) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
                else
                    RNGSNRDB(iSNR,iT,iNSC,iOf) = RNGSNRMIN;                 % ITER 1: Smallest SNR from range
                end
                ERR(1:NTX,iSNR,iT,iNSC,iOf) = 0;                            
                BITSTART = 1; IBITSTX = 0;
                while(BITSTART <= TOTALBITS-BPS(iNSC,iOf)+1)                % LOOP START MONTE CARLO
                    BITSTOP = BITSTART + BPS(iNSC,iOf) - 1;
                    sym_bits = randi([0 1],[NTX,BPS(iNSC,iOf)]);            % Generate random bits

                    % Generate OFDM Frame
                    for iTx = 1:NTX
                        sym_dec = bin2decMat(reshape(sym_bits(iTx,:),d,log2(M)))+1;
                        X(iTx,:) = genOFDMsignal(...                        
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
                    SNR = power(10,RNGSNRDB(iSNR,iT,iNSC,iOf)/10);
                    SNRrt = sqrt(SNR);
                    W = (XAVG(iT)/SNRrt)*randn(NRX,MODNSC);

                    % Compute output
                    Y = H*Xs + W;

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
                        ERR(iTx,iSNR,iT,iNSC,iOf) = ERR(iTx,iSNR,iT,iNSC,iOf) + errs;
                    end % Transmitters

                    IBITSTX = IBITSTX + BPS(iNSC,iOf);
                    BITSTART = BITSTOP + 1;
                end                                                         % LOOP STOP MONTE CARLO
                if (iT==2) && (iNSC==1) && (iOf==2)
                    abcd = 1;
                end
                BER(:,iSNR,iT,iNSC,iOf) = ERR(:,iSNR,iT,iNSC,iOf)/IBITSTX;
                    
                % Calcculate change in BER and SNR
                if iSNR>1
                    DBER(:,iSNR,iT,iNSC,iOf) = BER(:,iSNR,iT,iNSC,iOf) - BER(:,iSNR-1,iT,iNSC,iOf);
                    DSNR(iSNR,iT,iNSC,iOf) = RNGSNRDB(iSNR,iT,iNSC,iOf) - RNGSNRDB(iSNR-1,iT,iNSC,iOf);
                else
                    DBER(:,iSNR,iT,iNSC,iOf) = BER(:,iSNR,iT,iNSC,iOf);
                    DSNR(iSNR,iT,iNSC,iOf) = 0;
                end
                
                % Update progress on wait bar
                LOOPCOUNT = LOOPCOUNT + DSNR(iSNR,iT,iNSC,iOf);
                PROGRESS = LOOPCOUNT/TOTALLOOPS;
                waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...',PROGRESS*100));
                if(getappdata(hWB,'canceling'))
                    delete(hWB);
                    error('Simulation aborted');
                end
                
                % check and break if BER for ALL channels are below set thresholds
                if (RNGSNRDB(iSNR,iT,iNSC,iOf) > RNGSNRMAX) || isequal(sum(BER(:,iSNR,iT,iNSC,iOf) < BERTHMIN),NTX)
                    LOOPDONE = true;
                    LOOPCOUNT = (iT-1)*LENMODNSC*LENOFDMTYPES*RNGSNRLOOP +...
                                (iNSC-1)*LENOFDMTYPES*RNGSNRLOOP +...
                                (iOf-1)*RNGSNRLOOP + RNGSNRLOOP;
                    PROGRESS = LOOPCOUNT/TOTALLOOPS;
                    waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...',PROGRESS*100));
                    if(getappdata(hWB,'canceling'))
                        delete(hWB);
                        error('Simulation aborted');
                    end
                end
                iSNR = iSNR + 1;
            end % SNR                                                       % LOOP STOP SNR
        end % OFDM                                                          % LOOP STOP OFDM
    end % NSC                                                               % LOOP STOP NSC
end % CCT                                                                   % LOOP STOP CCT

if fSAVEALL                                                                 % SAVE script and data
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end
delete(hWB);
catch ex
delete(hWB);
rethrow(ex);
end

scrOFDMWDMPL(ctFileVars);                                                   % Call Show Results script





































