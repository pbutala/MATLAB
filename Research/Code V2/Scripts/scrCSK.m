% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;
addpath(genpath('..'));

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = true;
fCLOSEALL = false;
fARCHIVE = false;

fDECODER = 3; % 1.XYZ 2.RGB 3.TRIs
fCBC = 9; % 1<=fCBC<=9

fSHOWPGBAR = isequal(strfind(pwd,'graduate/pbutala'),[]);
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = '0_Scratch_';
% STRPREFIX = '1_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
csk = cCSK(fCBC);
RMN = csk.CBC(1).Center; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Band-i (~Red)
GMN = csk.CBC(2).Center; GSC = 1; GSD = 0;      % Mean, SD and scale to generate SPD of Band-j (~Green)
BMN = csk.CBC(3).Center; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Band-k (~Blue)
% switch fCBC
%     case 1
%         GMN = 564; GSC = 1; GSD = 0;      % Mean, SD and scale to generate SPD of Band-j (~Green)
%     case 2
%         GMN = 509; GSC = 1; GSD = 0;      % Mean, SD and scale to generate SPD of Band-j (~Green)
%     otherwise
%         error('CBC-%d not supported.',fCBC);
% end

cieFile = 'CIE1931_JV_1978_2deg';                 % CIE XYZ CMF curves file
flCIE = [cieFile '.csv'];
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
SPDTYP = SPDTYPE.GAUSSIAN;
cResp = cResponsivity();                    % Responsivity class instance
RResp = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
GResp = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
BResp = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)

NTX = 3; NRX = 3;

TOTALBITS = 1e4;                            % Total bit for transmtter to simulate

WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running CSK Simulation...'; % Wait Box title
FIGTITLE = 'Off';
MKTYP = {'o','+','^','s'};
% MKCLR = {'g','y','b','r'};
MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};

%% M-CSK CONSTELLATION
M = 4;
%      00,    01,    10,    11     (g w b r)
x = [csk.CBC(2).x 0 csk.CBC(3).x csk.CBC(1).x];
y = [csk.CBC(2).y 0 csk.CBC(3).y csk.CBC(1).y];
x(2) = mean(x([1 3 4]));
y(2) = mean(y([1 3 4]));
% switch fCBC
%     case 1
%         x = [0.402; 0.435; 0.169; 0.734];
%         y = [0.597; 0.290; 0.007; 0.265];
%     case 2
%         x = [0.011; 0.305; 0.169; 0.734];
%         y = [0.733; 0.335; 0.007; 0.265];
%     otherwise
%         error('CBC-%d not supported.',fCBC);
% end
Yc = 1;

%% ranges
RNGSNRMIN = 0; RNGSNRMAX = 100; SNROFST = 0;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000]; % DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR
BERTH = 1e-3;   BERTHMIN = 0.75*BERTH;       % BER thresholds;

RNGSNRST = 20:10:RNGSNRMAX;   LENSNRST = numel(RNGSNRST);
RNGBERST = power(10,[-3 -4]);     LENBERST = numel(RNGBERST);
IDXSNRST = 1; IDXBERST = 1; IDXCHST = 0;
%% config

ctDirRes = '..\..\..\..\MatlabResults\14. CSK\';
ctDirData = [ctDirRes STRPREFIX 'Data\'];
ctFileCodeSrc = '.\scrCSK.m';
ctFileCodeDest = [ctDirData STRPREFIX 'scrCSK' CHARIDXARCHIVE '.m']; % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];   % Data file name
ctFileChnlStPRE = [ctDirData STRPREFIX 'datChnlStat'];   % Channel state file name
if ~exist(ctDirData,'dir')   % if data directory does NOT exist
    mkdir(ctDirData);        % create data dir
end

%% logic
try
    delete([ctDirData '*' CHAROVERWRITE '.*']);
    %% prep start
    % Wait Bar to show progress
    if fSHOWPGBAR
        hWB = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
            'Visible','Off',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
        setappdata(hWB,'canceling',0);
    end
    LOOPCOUNT = 1; 
    TOTALLOOPS = RNGSNRLOOP;
    
    %% compute
    switch SPDTYP
        case SPDTYPE.GAUSSIAN 
            sSPDTYP = 'Gaussian';
        case SPDTYPE.LORENTZIAN
            sSPDTYP = 'Lorentzian';
        otherwise
            error('SPDTYPE must be either ''Gaussian'' or ''Lorentzian''');
    end
    ctMatDir = '..\Matfiles\LEDPSD\';
    sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
        RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
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
        RGB.initialize(fSHOWPGBAR);                   % Initialize led
        if ~exist(sPSDDIR,'dir')
            mkdir(sPSDDIR);                 % Create directory to store characterization data
        end
        save(RGBledmat,'RGB');              % Save LED characterization
    end    
    
    BITERR = 0;
    RNGSNRDB = [];
    
    %% Generate Symbols
    SYMS = zeros(NTX,M);
    TRIS = zeros(NTX,M);
    PTXAVG = 0;
    for iSym = 1:M
        xc = x(iSym); yc = y(iSym);
        [S,Ds,Ts] = RGB.getPSD(xc,yc);       % Get transmitter(s) SPD for set x,y
%         SYMS(:,iSym) = Ts;
        for iTx = 1:NTX
            TRIS(iTx,iSym) = Ts(iTx);
            SYMS(iTx,iSym) = Ds{iTx}.rdFlux;
        end
        PTXAVG = PTXAVG + S.rdFlux/M;
    end
    
    %% Compute channel matrix
    H = [RResp 0 0;0 GResp 0; 0 0 BResp];
    
    %% Compute SYMS at Receiver
    SYMSRX = zeros(NRX,M);
    for iSym = 1:M
        SYMSRX(:,iSym) = H*SYMS(:,iSym);
    end
    
    TSTART = tic;
    BITSSYM = log2(M);
    LOOPDONE = false; iSNR = 1;
    while ~LOOPDONE                                                 % LOOP START SNR (dynamically select next SNR)
        if iSNR > 1
            dber = BER(iSNR-1);   % Smallest BER from all those above the BERTH limit
            % dynamically select next SNR based on how close the current BER is to threshold
            RNGSNRDB(iSNR) = RNGSNRDB(iSNR-1) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
        else
            RNGSNRDB(iSNR) = RNGSNRMIN;                 % ITER 1: Smallest SNR from range
        end
        
        SNR = power(10,RNGSNRDB(iSNR)/10);
        SNRrt = sqrt(SNR);
        BITERR = 0; BITCOUNT = 0;
        ARLEN = floor(TOTALBITS/BITSSYM);
        SZRXSYMS = [3 ARLEN];
        switch fDECODER
            case 1
                SYMSEST = SYMS;
            case 2
                SYMSEST = [x,y,1-x-y]';
            case 3
                SYMSEST = TRIS;
        end
        CHST = cChnlState(ARLEN,SZRXSYMS,SYMSEST,H);
        CHST.SNRdB = RNGSNRDB(iSNR);
        MCIDX = 0;
        while(BITCOUNT < TOTALBITS - BITSSYM + 1)
            MCIDX  = MCIDX + 1;                         % increment monte-carlo index
            BITS = randi([0 1],[1 BITSSYM]);            % Generate random bits
            SYMIDX = bin2decMat(BITS)+1;
            CHST.TxIdx(MCIDX) = SYMIDX;
            
            X = SYMS(:,SYMIDX);
            
            % Compute noise
            W = (PTXAVG/SNRrt)*randn(NRX,1);
            Y = H*X + W;
            
            %% Estimate transmitted vector
            Xh = H\Y;
            
            switch fDECODER
                case 1
                    % MAP estimate (on RGB data) TODO: Has Errors. Need to Fix
                    for iSym = 1:M
                        vYSym = Y-SYMSRX(:,iSym);
                        dYSym(iSym) = sum((vYSym.*vYSym),1);
                    end
                    CHST.RxSymEst(:,MCIDX) = Xh;
                case 2
                    % MAP estimate (on xyz data) TODO: Has Errors. Need to Fix
                    XYZh = RGB.Txyz*Xh;
                    xyzh = XYZh/sum(XYZh,1);
                    vYSym = zeros(3,1);
                    for iSym = 1:M
                        vYSym(1) = xyzh(1)-x(iSym);
                        vYSym(2) = xyzh(2)-y(iSym);
                        dYSym(iSym) = sum((vYSym.*vYSym),1);
                    end
                    CHST.RxSymEst(:,MCIDX) = xyzh;
                case 3
                    % MAP estimate (on tristimulus values)
                    % calibrate for primary powers
                    th(1) = Xh(1)/RGB.PSDs{1}.rdFlux;
                    th(2) = Xh(2)/RGB.PSDs{2}.rdFlux;
                    th(3) = Xh(3)/RGB.PSDs{3}.rdFlux;
                    th = th/sum(th);
                    
                    vYSym = zeros(3,1);
                    for iSym = 1:M
                        vYSym(1) = th(1)-TRIS(1,iSym);
                        vYSym(2) = th(2)-TRIS(2,iSym);
                        vYSym(3) = th(3)-TRIS(3,iSym);
                        dYSym(iSym) = sum((vYSym.*vYSym),1);
                    end
                    CHST.RxSymEst(:,MCIDX) = th;
                otherwise
                    error('Decoder not defined');
            end
%             SYMIDXh = find(dYSym == min(dYSym),1,'first');
            [~,SYMIDXh] = min(dYSym,[],2);
            CHST.RxIdx(MCIDX) = SYMIDXh;
            
            BITSh = dec2binMat(SYMIDXh-1,log2(M));

            BITERR = BITERR + biterr2(BITS, BITSh);
            BITCOUNT = BITCOUNT + BITSSYM;
        end
        
        BER(iSNR) = BITERR/BITCOUNT;
        CHST.BER = BER(iSNR);
        
        % Determine if state should be saved or not
        FLGST = false;
        if (IDXSNRST <= LENSNRST)
            if (RNGSNRDB(iSNR) > RNGSNRST(IDXSNRST))
                FLGST = true;
                IDXSNRST = IDXSNRST + 1;
            end
        end
        
        if (IDXBERST <= LENBERST)
            if (BER(iSNR) < RNGBERST(IDXBERST))
                FLGST = true;
                IDXBERST = IDXBERST + 1;
            end
        end
        
        if FLGST
            IDXCHST = IDXCHST + 1;
            FileChnlSt = [ctFileChnlStPRE sprintf('%d',IDXCHST) CHARIDXARCHIVE '.mat'];
            save(FileChnlSt,'CHST');
        end
        clear CHST;
        % Calculate change in BER and SNR
        if iSNR>1
            DSNR = RNGSNRDB(iSNR) - RNGSNRDB(iSNR-1);
        else
            DSNR = 0;
        end
        % Update progress on wait bar
        LOOPCOUNT = LOOPCOUNT + DSNR;
        PROGRESS = LOOPCOUNT/TOTALLOOPS;
        TELAPSED = toc(TSTART);
        TREM = (TELAPSED/PROGRESS)-TELAPSED;
        if fSHOWPGBAR
            waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
            if(getappdata(hWB,'canceling'))
                delete(hWB);
                error('Simulation aborted');
            end
        end
                            
        % check and break if BER for ALL channels are below set thresholds
        if (RNGSNRDB(iSNR) >= RNGSNRMAX) || (BER(iSNR) < BERTHMIN)
            LOOPDONE = true;
            LOOPCOUNT = RNGSNRLOOP;
            PROGRESS = LOOPCOUNT/TOTALLOOPS;
            TELAPSED = toc(TSTART);
            TREM = (TELAPSED/PROGRESS)-TELAPSED;
            if fSHOWPGBAR
                waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %s',PROGRESS*100,getTimeString(TREM)));
                if(getappdata(hWB,'canceling'))
                    delete(hWB);
                    error('Simulation aborted');
                end
            end
        end
        iSNR = iSNR + 1;
    end
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    %% prep stop
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
catch ex
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    %     setpref('Internet','E_mail','pbutala@bu.edu');
    %     setpref('Internet','SMTP_Server','smtp.bu.edu');
    %     STREMAIL = ['Simulation ' STRPREFIX ' error!'];
    %     sendmail('pankil.butala@gmail.com',STREMAIL);
    rethrow(ex);
end
% setpref('Internet','E_mail','pbutala@bu.edu');
% setpref('Internet','SMTP_Server','smtp.bu.edu');
% STREMAIL = ['Simulation ' STRPREFIX ' done. Starting plots.'];
% sendmail('pankil.butala@gmail.com',STREMAIL);
fprintf('--scrCSK Done--\n');
scrCSKPL;                                                  % Call Show Results script






























































