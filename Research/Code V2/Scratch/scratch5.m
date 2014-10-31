% scrCSKOFDM
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% FLAGS
% fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = true;
fCLOSEALL = true;
fARCHIVE = true;
fDECODER = 3; % 1.XYZ 2.RGB 3.TRIs
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = 'Scratch_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

% CIE1931 RGB monochromatic wavelengths
% RMN = 700; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Red led
% GMN = 546; GSC = 1; GSD = 0;               % Mean, SD and scale to generate SPD of Green led
% BMN = 436; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led

RMN = 703; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Red led
GMN = 564; GSC = 1; GSD = 0;               % Mean, SD and scale to generate SPD of Green led
% BMN = 429; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led
BMN = 509; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led
cieFile = 'CIE1931_JV_1978_2deg';                 % CIE XYZ CMF curves file
flCIE = [cieFile '.csv'];
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
SPDTYP = SPDTYPE.GAUSSIAN;
cResp = cResponsivity();                    % Responsivity class instance
RResp = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
GResp = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
BResp = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)

NTX = 3; NRX = 3;

TOTALBITS = 1e3;                            % Total bit for transmtter to simulate

WBTITLE = 'Running CSK-OFDM Simulation...'; % Wait Box title
FIGTITLE = 'Off';
MKTYP = {'o','+','^','s'};
MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};

%% M-CSK CONSTELLATION
C_M = 4;
%      00,    01,    10,    11     (g w b r)
% CBC 1
% x = [0.402; 0.435; 0.169; 0.734];
% y = [0.597; 0.290; 0.007; 0.265];
% CBC 2: this is as specified in standard. varies slightly from actual
% calculation. !!                                                           **********IMP**********
x = [0.402; 0.382; 0.011; 0.734];
y = [0.597; 0.532; 0.733; 0.265];
Yc = 1;

%% ranges
RNGSNRMIN = 0; RNGSNRMAX = 150; SNROFST = 0;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000]; 
% DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR
BERTH = 1e-2;   BERTHMIN = 0.5*BERTH;       % BER thresholds;

% OFDM RANGES
RNGOFDMTYPES = {'dcoofdm';'acoofdm'};   LENOFDMTYPES = numel(RNGOFDMTYPES); % OFDM types
% RNGOFDMOFST = [3.2 3.2];                  LENOFDMOFST = size(RNGOFDMOFST,2);  % OFDM offsets
RNGMODDCO = power(2,2);               LENMOD  = numel(RNGMODDCO);           % Subcarrier modulation order
RNGMODACO = RNGMODDCO.^2;
RNGMODNSC = power(2,6);                 LENMODNSC = numel(RNGMODNSC);       % Number of subcarriers

DOFST = 2;
RNGOFDMOFSTACOXTR = 0.2;
RNGOFDMOFSTDCOXTR = 3.2;                              LENOFSTIGNR = numel(RNGOFDMOFSTDCOXTR);   % OFDM extra offsets  
RNGOFDMOFSTACO = [0:DOFST:4 RNGOFDMOFSTACOXTR];       
RNGOFDMOFSTDCO = [0:DOFST:4 RNGOFDMOFSTDCOXTR];       LENOFDMOFST = numel(RNGOFDMOFSTDCO);      % OFDM offsets

% RNGOFDMOFSTACO = 0.2;
% RNGOFDMOFSTDCO = 3.2;         LENOFDMOFST = numel(RNGOFDMOFSTDCO);  % OFDM offsets

% O_RES = 1e-3; O_PERS = power(2,10); O_MAXMC = 1e9;
WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
%% config

% STATION
ctDirRes = '..\..\..\..\MatlabResults\15. CSKOFDM\';
ctDirData = [ctDirRes STRPREFIX 'Data\'];
ctFileCodeSrc = '.\scratch2.m';
        
% switch fSTATION
%     % Results directory; Spource file; LED table dir;
%     case 1
%         ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\15. CSKOFDM\';
%         ctDirData = [ctDirRes STRPREFIX 'Data\'];
%         ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrCSKOFDM.m';
%     case 2
%         ctDirRes = '/home/pbutala/My Documents/MatlabResults/15. CSKOFDM/';
%         ctDirData = [ctDirRes STRPREFIX 'Data/'];
%         ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrCSKOFDM.m';
%     case 3
%         ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\15. CSKOFDM\\';
%         ctDirData = [ctDirRes STRPREFIX 'Data\\'];
%         ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSKOFDM.m';
%     case 4
%         ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\15. CSKOFDM\\';
%         ctDirData = [ctDirRes STRPREFIX 'Data\\'];
%         ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSKOFDM.m';
%     otherwise
%         error('Station not defined');
% end
ctFileCodeDest = [ctDirData STRPREFIX 'scrCSKOFDM' CHARIDXARCHIVE '.m']; % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datCSKOFDM' CHARIDXARCHIVE '.mat'];   % Data file name
if ~exist(ctDirData,'dir')   % if data directory does NOT exist
    mkdir(ctDirData);        % create data dir
end

%% logic
try
    delete([ctDirData '*' CHAROVERWRITE '.*']);
    %% prep start
    % Wait Bar to show progress
    hWB = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
        'Visible','Off',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
    setappdata(hWB,'canceling',0);
    LOOPCOUNT = 1;
    TOTALLOOPS = RNGSNRLOOP*LENMOD*LENMODNSC*LENOFDMOFST*LENOFDMTYPES;
    
    %% compute
    switch SPDTYP
        case SPDTYPE.GAUSSIAN
            sSPDTYP = 'Gaussian';
        case SPDTYPE.LORENTZIAN
            sSPDTYP = 'Lorentzian';
        otherwise
            error('SPDTYPE must be either ''Gaussian'' or ''Lorentzian''');
    end
    % STATION
    ctMatDir = '..\Matfiles\LEDPSD\';
    sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
        RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
            
%     switch fSTATION
%         % Results directory; Spource file; LED table dir;
%         case 1
%             ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\';
%             sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
%                 RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
%         case 2
%             ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/LEDPSD/';
%             sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '/' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
%                 RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '/'];
%         case 3
%             ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
%             sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
%                 RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
%         case 4
%             ctMatDir = 'C:\\Users\\Pankil\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
%             sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
%                 RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
%         otherwise
%             error('Station not defined');
%     end
    RGBledmat = [sPSDDIR sprintf('res_%0.5f',RES) '.mat'];                  % LED table mat-file
    %% SPDs
    switch SPDTYP
        case SPDTYPE.GAUSSIAN
            Rspd = getSOG(RMN,RSD,RSC,lambdas);
            Gspd = getSOG(GMN,GSD,GSC,lambdas);
            Bspd = getSOG(BMN,BSD,BSC,lambdas);
        case SPDTYPE.LORENTZIAN
            Rspd = getLorentzian(RMN,RSD,lambdas);
            Gspd = getLorentzian(GMN,GSD,lambdas);
            Bspd = getLorentzian(BMN,BSD,lambdas);
        otherwise
            error('SPDTYPE must be either ''Gaussian'' or ''Lorentzian''');
    end
    
    Rspd = Rspd/(sum(Rspd)*LAMBDADELTA);
    Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rspd,'Red');   % Normalized RED SPD
    Gspd = Gspd/(sum(Gspd)*LAMBDADELTA);
    Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gspd,'Green');   % Normalized GREEN SPD
    Bspd = Bspd/(sum(Bspd)*LAMBDADELTA);
    Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bspd,'Blue');   % Normalized BLUE SPD
    
    %% variables
    if exist(RGBledmat,'file')              % If LED characterization table exists
        load(RGBledmat,'RGB');              % Load RGB led
    else
        RGB = cLEDs(RES,Rch,Gch,Bch,flCIE);
        RGB.initialize();                   % Initialize led
        if ~exist(sPSDDIR,'dir')
            mkdir(sPSDDIR);                 % Create directory to store characterization data
        end
        save(RGBledmat,'RGB');              % Save LED characterization
    end
    
    RNGSNRDB = [];
    
    %% Generate CSK Symbols
    C_SYMS = zeros(NTX,C_M);
    C_TRIS = zeros(NTX,C_M);
    PTXAVG = 0; PTXAVGh = 0;
    for iSym = 1:C_M
        xc = x(iSym); yc = y(iSym);
        [S,Ds,Ts] = RGB.getPSD(xc,yc);       % Get transmitter(s) SPD for set x,y
%         C_SYMS(:,iSym) = Ts;
        for iTx = 1:NTX
            C_TRIS(iTx,iSym) = Ts(iTx);
            C_SYMS(iTx,iSym) = Ds{iTx}.rdFlux;
        end
        PTXAVG = PTXAVG + S.rdFlux/C_M;
    end
    
    %% Compute channel matrix
    H = [RResp 0 0;0 GResp 0; 0 0 BResp];
    
    %% Compute SYMS at Receiver
    C_SYMSRX = zeros(NRX,C_M);
    for iSym = 1:C_M
        C_SYMSRX(:,iSym) = H*C_SYMS(:,iSym);
    end
    
    TSTART = tic;
    LOOPIDX = 0;
    % Loop through different configurations
    for iM = 1:LENMOD
        for iNSC = 1:LENMODNSC                                                  % LOOP START MODNSC
            MODNSC = RNGMODNSC(iNSC);
            C_BPS(iNSC) = MODNSC*log2(C_M);
            for iDC = 1:LENOFDMOFST
                for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
                    LOOPIDX = LOOPIDX + 1;
                    ofdmType = lower(RNGOFDMTYPES{iOf});
                    switch lower(ofdmType)
                        case 'acoofdm'
                            d = MODNSC/4;
                            STROFDM = 'ACO';
                            OffsetDcoStddev = 0;
                            OffsetAcoStddev = RNGOFDMOFSTACO(iDC);
                            O_M = RNGMODACO(iM);
                        case 'dcoofdm'
                            d = MODNSC/2-1;
                            STROFDM = 'DCO';
                            OffsetDcoStddev = RNGOFDMOFSTDCO(iDC);
                            OffsetAcoStddev = 0;
                            O_M = RNGMODDCO(iM);
                    end
                    O_SYMS = getQAMsyms(O_M);                                           % Get QAM symbols
                    O_BPS(iM,iNSC,iOf) = d*log2(O_M);                                      % Calculate Bits Per Symbol
                    
                    waitbar(0,hWB,sprintf('Computing mean OFDM signal...'));
                    [~,~,O_tSig_AVG(iM,iNSC,iDC,iOf)] = getOFDMMean(ofdmType, MODNSC, O_M, OffsetDcoStddev, OffsetAcoStddev);
                    LOOPDONE = false; iSNR = 1;
                    while ~LOOPDONE                                                 % LOOP START SNR (dynamically select next SNR)
                        if iSNR > 1
                            dber = BER(iSNR-1,iM,iNSC,iDC,iOf);   % Smallest BER from all those above the BERTH limit
                            % dynamically select next SNR based on how close the current BER is to threshold
                            RNGSNRDB(iSNR,iM,iNSC,iDC,iOf) = RNGSNRDB(iSNR-1,iM,iNSC,iDC,iOf) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
                        else
                            RNGSNRDB(iSNR,iM,iNSC,iDC,iOf) = RNGSNRMIN;                    % ITER 1: Smallest SNR from range
                        end
                        
                        SNR = power(10,RNGSNRDB(iSNR,iM,iNSC,iDC,iOf)/10);
                        SNRrt = sqrt(SNR);
                        BITERR = 0; C_BITERR = 0; O_BITERR = 0;
                        BITCOUNT = 0; C_BITCOUNT = 0; O_BITCOUNT = 0;
                        
                        MCIDX = 0; 
                        PTXAVGh(iSNR,iM,iNSC,iDC,iOf) = 0;
                        O_tSig_AVGh(iSNR,iM,iNSC,iDC,iOf) = 0;
                        while(BITCOUNT < TOTALBITS - C_BPS(iNSC) - O_BPS(iM,iNSC,iOf) + 1)
                            MCIDX  = MCIDX + 1;                                     % increment monte-carlo index
                            
                            C_BITS = randi([0 1],[MODNSC log2(C_M)]);                        % Generate CSK bits
                            C_SYMIDX = bin2decMat(C_BITS)+1;
                            
                            O_BITS = randi([0 1],[d log2(O_M)]);              % Generate OFDM bits
                            O_SYMIDX = bin2decMat(O_BITS)+1;
                            
                            tSig = genOFDMsignal(...
                                'data',O_SYMIDX,...
                                'OFDMtype',ofdmType,...
                                'N',MODNSC,...
                                'Symbols',O_SYMS,...
                                'OffsetDcoStddev', OffsetDcoStddev,...
                                'OffsetAcoStddev', OffsetAcoStddev)';
                            tSigMat = repmat(tSig,NTX,1);
                            O_tSig_AVGh(iSNR,iM,iNSC,iDC,iOf) = (O_tSig_AVGh(iSNR,iM,iNSC,iDC,iOf)*(MCIDX-1) + mean(tSig))/MCIDX;
                            
                            X = C_SYMS(:,C_SYMIDX).*tSigMat;
                            X = X/O_tSig_AVG(iM,iNSC,iDC,iOf);                         % Normalize signal to unit Mean (average) signal.
                            
                            PTXAVGh(iSNR,iM,iNSC,iDC,iOf) = (PTXAVGh(iSNR,iM,iNSC,iDC,iOf)*(MCIDX-1) + sum(X(:))/MODNSC)/MCIDX;
                            % Compute noise
                            W = (PTXAVG/SNRrt)*randn(NRX,MODNSC);
                            Y = H*X + W;
                            % Y = H*X;
                            
                            %% Estimate transmitted vector
                            Xh = H\Y;
                            
%                             %% calibrate for primary powers
%                             Xh(1,:) = Xh(1,:)/RGB.PSDs{1}.rdFlux;
%                             Xh(2,:) = Xh(2,:)/RGB.PSDs{2}.rdFlux;
%                             Xh(3,:) = Xh(3,:)/RGB.PSDs{3}.rdFlux;
%                             
%                             %% Estimate CSK frame
%                             XYZh = RGB.Txyz*Xh;
%                             xyzh = XYZh./repmat(sum(XYZh,1),3,1);
                            
                            switch fDECODER
                                case 1
                                    % MAP estimate (on RGB data) TODO: Has Errors. Need to Fix
                                    for iSym = 1:C_M
                                        vYSym = Y-C_SYMSRX(:,iSym);
                                        dYSym(iSym) = sum((vYSym.*vYSym),1);
                                    end
                                case 2
                                    %% MAP estimate (on xyz data) TODO: Has Errors. Need to Fix
                                    vYSym = zeros(2,MODNSC);
                                    for iSym = 1:C_M
                                        vYSym(1,:) = xyzh(1,:)-x(iSym);
                                        vYSym(2,:) = xyzh(2,:)-y(iSym);
                                        dYSym(iSym,:) = sum((vYSym.*vYSym),1);
                                    end
                                case 3
                                    %% MAP estimate (on tristimulus values)
                                    % calibrate for primary powers
                                    th(1,:) = Xh(1,:)/RGB.PSDs{1}.rdFlux;
                                    th(2,:) = Xh(2,:)/RGB.PSDs{2}.rdFlux;
                                    th(3,:) = Xh(3,:)/RGB.PSDs{3}.rdFlux;
                                    th = th./repmat(sum(th,1),3,1);
                                    
                                    vYSym = zeros(3,MODNSC);
                                    for iSym = 1:C_M
                                        vYSym(1,:) = th(1,:)-C_TRIS(1,iSym);
                                        vYSym(2,:) = th(2,:)-C_TRIS(2,iSym);
                                        vYSym(3,:) = th(3,:)-C_TRIS(3,iSym);
                                        dYSym(iSym,:) = sum((vYSym.*vYSym),1);
                                    end
                                otherwise
                                    error('Decoder not defined');
                                    
                            end
                            [~,C_SYMIDXh] = min(dYSym,[],1);
                            C_BITSh = dec2binMat(C_SYMIDXh-1,log2(C_M));
                            C_BITERRd = biterr2(C_BITS, C_BITSh);
                            C_BITERR = C_BITERR + C_BITERRd;
                            C_BITCOUNT = C_BITCOUNT + C_BPS(iNSC);
                            
                            %% Estimate OFDM frame
                            tSigh = zeros(1,MODNSC);
                            for iSPL = 1:MODNSC
                                tSigh(iSPL) = C_SYMS(:,C_SYMIDXh(iSPL))\Xh(:,iSPL);
                            end
                            
                            %% Scale tSigh based on average signal value
                            tSigh = tSigh*O_tSig_AVG(iM,iNSC,iDC,iOf); 
                            
                            O_SYMIDXh = decodeOFDMsignal(tSigh,...
                                'OFDMtype',ofdmType,...
                                'N',MODNSC,...
                                'Symbols',O_SYMS);
                            
                            O_BITSh = dec2binMat(O_SYMIDXh-1,log2(O_M));
                            O_BITERRd = biterr2(O_BITS, O_BITSh);
                            O_BITERR = O_BITERR + O_BITERRd;
                            O_BITCOUNT = O_BITCOUNT + O_BPS(iM,iNSC,iOf);
                            
                            BITERR = C_BITERR + O_BITERR;
                            BITCOUNT = C_BITCOUNT + O_BITCOUNT;
                            
%                             if (C_BITERRd > 0) || (O_BITERRd > 0)
%                                 fprintf('Error in %s frame\n',ofdmType);
%                             end
                        end
                        
                        C_BER(iSNR,iM,iNSC,iDC,iOf) = C_BITERR/C_BITCOUNT;
                        O_BER(iSNR,iM,iNSC,iDC,iOf) = O_BITERR/O_BITCOUNT;
                        BER(iSNR,iM,iNSC,iDC,iOf) = BITERR/BITCOUNT;
                        
                        % Calculate change in BER and SNR
                        if iSNR>1
                            DSNR = RNGSNRDB(iSNR,iM,iNSC,iDC,iOf) - RNGSNRDB(iSNR-1,iM,iNSC,iDC,iOf);
                        else
                            DSNR = 0;
                        end
                        % Update progress on wait bar
                        LOOPCOUNT = LOOPCOUNT + DSNR;
                        PROGRESS = LOOPCOUNT/TOTALLOOPS;
                        TELAPSED = toc(TSTART);
                        TREM = (TELAPSED/PROGRESS)-TELAPSED;
                        waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
                        if(getappdata(hWB,'canceling'))
                            delete(hWB);
                            error('Simulation aborted');
                        end
                        
                        % check and break if BER for ALL channels are below set thresholds
                        if (RNGSNRDB(iSNR,iM,iNSC,iDC,iOf) >= RNGSNRMAX) || (BER(iSNR,iM,iNSC,iDC,iOf) < BERTHMIN)
                            LOOPDONE = true;
                            % LOOPCOUNT = RNGSNRLOOP*iM*iNSC*iDC*iOf;
                            LOOPCOUNT = RNGSNRLOOP*LOOPIDX;
                            PROGRESS = LOOPCOUNT/TOTALLOOPS;
                            TELAPSED = toc(TSTART);
                            TREM = (TELAPSED/PROGRESS)-TELAPSED;
                            waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
                            if(getappdata(hWB,'canceling'))
                                delete(hWB);
                                error('Simulation aborted');
                            end
                        end
                        iSNR = iSNR + 1;
                    end                                                         % WHILE
                end                                                             % LOOP STOP OFDM types
            end                                                                 % LOOP STOP OFFSET
        end                                                                     % LOOP STOP MODNSC
    end
    fprintf('Loops done\n');
    save(ctFileVars);                       % save workspace
    fprintf('Saved workspace (init)\n');
    I = find(RNGSNRDB == 0);
%     [D1,D2,D3,D4,D5] = ind2sub(size(RNGSNRDB),I);
%     RNGSNRDB(D1(D1~=1),D2(D1~=1),D3(D1~=1),D4(D1~=1),D5(D1~=1)) = nan;
%     BER(D1(D1~=1),D2(D1~=1),D3(D1~=1),D4(D1~=1),D5(D1~=1)) = nan;
    RNGSNRDB(I) = nan;
    BER(I) = nan;
    %% prep stop
    fprintf('Saving script and workspace\n');
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
    delete(hWB);
catch ex
    delete(hWB);
%         setpref('Internet','E_mail','pbutala@bu.edu');
%         setpref('Internet','SMTP_Server','smtp.bu.edu');
%         STREMAIL = ['Simulation ' STRPREFIX ' error!'];
%         sendmail('pankil.butala@gmail.com',STREMAIL);
    rethrow(ex);
end
% setpref('Internet','E_mail','pbutala@bu.edu');
% setpref('Internet','SMTP_Server','smtp.bu.edu');
% STREMAIL = ['Simulation ' STRPREFIX ' done. Starting plots.'];
% sendmail('pankil.butala@gmail.com',STREMAIL);
% fprintf('--scrCSKOFDM Done--\n');
% scrCSKOFDMPL;                                                  % Call Show Results script






























































