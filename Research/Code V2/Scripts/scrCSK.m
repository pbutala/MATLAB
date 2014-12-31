% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% config
RNGCBC = 1:9;                               % CBCs to consider
M = 16;
TOTALBITS = 1e4;                            % Total bit for transmtter to simulate
% DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR

% FLAGS
fSAVEALL = true;
fCLOSEALL = true;
fSAVECHST = false;
fDECODER = 2; % 1.RGB 2.XYZ 3.TRIs
fSHOWPGBAR = isequal(strfind(pwd,'graduate/pbutala'),[]);
fARCHIVE = true;
CHAROVERWRITE = '~';
STRPREFIX = 'M16_1e4_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
rng('default');

%% set paths
ctFileCodeSrc = [mfilename('fullpath') '.m'];                           % get fullpath of current file
[ctScrDir,~,~] = fileparts(ctFileCodeSrc);                              % get scripts dir
cd(ctScrDir);                                                           % set scripts dir as pwd (reference)
ctDirRes = '..\..\..\..\MatlabResults\14. CSK\';
ctDirData = [ctDirRes STRPREFIX 'Data\'];
ctFileCodeDest = [ctDirData STRPREFIX 'scrCSK' CHARIDXARCHIVE '.m'];    % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];      % Data file name
ctFileChnlStPRE = [ctDirData STRPREFIX 'datChnlStat'];                  % Channel state file name
if ~exist(ctDirData,'dir')                                              % if data directory does NOT exist
    mkdir(ctDirData);                                                   % create data dir
end
ctMatDir = '..\Matfiles\LEDPSD\';

addpath(genpath('..'));

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
TXLM = 1;                                           % Lumens per symbol to transmit
cieFile = 'CIE1931_JV_1978_2deg';                 % CIE XYZ CMF curves file
flCIE = [cieFile '.csv'];
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
SPDTYP = SPDTYPE.GAUSSIAN;
%% compute
switch SPDTYP
    case SPDTYPE.GAUSSIAN
        sSPDTYP = 'Gaussian';
    case SPDTYPE.LORENTZIAN
        sSPDTYP = 'Lorentzian';
    otherwise
        error('SPDTYPE must be either ''Gaussian'' or ''Lorentzian''');
end
cResp = cResponsivity();                    % Responsivity class instance

NTX = 3; NRX = 3;
NRXrt = sqrt(NRX);

WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running CSK Simulation...'; % Wait Box title
FIGTITLE = 'Off';
MKTYP = {'o','+','^','s'};
% MKCLR = {'g','y','b','r'};
MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};

%% ranges
LENCBC = numel(RNGCBC);                     % number of CBCs to consider

RNGSNRMIN = 0; RNGSNRMAX = 100; SNROFST = 0;
RNGSNRMINPL = 0; RNGSNRMAXPL = 60;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000];
BERTH = 1e-3;   BERTHMIN = 0.5*BERTH;       % BER thresholds;

if fSAVECHST
    RNGSNRST = 20:10:RNGSNRMAX;   LENSNRST = numel(RNGSNRST);
    RNGBERST = power(10,[-3 -4]);     LENBERST = numel(RNGBERST);
    IDXSNRST = ones(LENCBC,1);
    IDXBERST = ones(LENCBC,1);
end
IDXCHST = zeros(LENCBC,1); % count number of states saved per CBC. used in scrCSKPL
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
    TOTALLOOPS = LENCBC*RNGSNRLOOP;
    RNGSNRDB = [];
    PTXAVG = zeros(LENCBC,1);
    PAVGPERLM = zeros(LENCBC,1);
    % loop over CSK CBCs
    for iCBC = 1:LENCBC
        fCBC = RNGCBC(iCBC);
        
        % CSK object to get CBCs
        csk(iCBC) = cCSK(fCBC,M);                                       % CSK object
        RMN = csk(iCBC).CBC(1).Center; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Band-i (~Red)
        GMN = csk(iCBC).CBC(2).Center; GSC = 1; GSD = 0;      % Mean, SD and scale to generate SPD of Band-j (~Green)
        BMN = csk(iCBC).CBC(3).Center; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Band-k (~Blue)
        
        % Get responsivity for three bands
        RResp = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        GResp = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        BResp = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        
        [x,y] = csk(iCBC).getSyms();
        
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
        RGBLED(iCBC) = RGB;
        BITERR = 0;
        
        %% Generate Symbols
        SYMS = zeros(NTX,M);
        LFLX = zeros(1,M);
        TRIS = zeros(NTX,M);
        for iSym = 1:M
            xc = x(iSym); yc = y(iSym);
            [S,Ds,Ts] = RGBLED(iCBC).getPSD(xc,yc);       % Get transmitter(s) SPD for set x,y
            for iTx = 1:NTX
                TRIS(iTx,iSym) = Ts(iTx);
                SYMS(iTx,iSym) = Ds{iTx}.rdFlux;
                LFLX(iSym) = LFLX(iSym) + Ds{iTx}.lmFlux;
            end
            TRIS(:,iSym) = TXLM*TRIS(:,iSym)/LFLX(iSym);
            SYMS(:,iSym) = TXLM*SYMS(:,iSym)/LFLX(iSym);
            PTXAVG(iCBC) = PTXAVG(iCBC) + sum(SYMS(:,iSym),1)/M;
            PAVGPERLM(iCBC) = PTXAVG(iCBC)/TXLM;
        end
        XAVG = mean(SYMS,2);
        
        %% Compute channel matrix
        H = [RResp 0 0;0 GResp 0; 0 0 BResp];
        
        %% Compute average received signal power
        SIGRXAVG = sqrt(trace(H*XAVG*XAVG'*H'));
        
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
                dber = BER(iSNR-1,iCBC);   % Smallest BER from all those above the BERTH limit
                % dynamically select next SNR based on how close the current BER is to threshold
                RNGSNRDB(iSNR,iCBC) = RNGSNRDB(iSNR-1,iCBC) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
            else
                RNGSNRDB(iSNR,iCBC) = RNGSNRMIN;                 % ITER 1: Smallest SNR from range
            end
            
            SNR = power(10,RNGSNRDB(iSNR,iCBC)/10);
            SNRrt = sqrt(SNR);
            BITERR = 0; BITCOUNT = 0;
            if fSAVECHST
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
                CHST(iCBC) = cChnlState(ARLEN,SZRXSYMS,SYMSEST,H);
                CHST(iCBC).SNRdB = RNGSNRDB(iSNR,iCBC);
            end
            MCIDX = 0;
            while(BITCOUNT < TOTALBITS - BITSSYM + 1)
                MCIDX  = MCIDX + 1;                         % increment monte-carlo index
                BITS = randi([0 1],[1 BITSSYM]);            % Generate random bits
                SYMIDX = bin2decMat(BITS)+1;
                if fSAVECHST
                    CHST(iCBC).TxIdx(MCIDX) = SYMIDX;
                end
                X = SYMS(:,SYMIDX);
                
                % Compute noise
                Wsd = SIGRXAVG/SNRrt;
                W = (Wsd/NRXrt)*randn(NRX,1);
                
                % channel
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
                        if fSAVECHST
                            CHST(iCBC).RxSymEst(:,MCIDX) = Xh;
                        end
                    case 2
                        % MAP estimate (on xyz data)
                        XYZh = RGBLED(iCBC).Txyz*Xh;
                        xyzh = XYZh/sum(XYZh,1);
                        vYSym = zeros(3,1);
                        for iSym = 1:M
                            vYSym(1) = xyzh(1)-x(iSym);
                            vYSym(2) = xyzh(2)-y(iSym);
                            dYSym(iSym) = sum((vYSym.*vYSym),1);
                        end
                        if fSAVECHST
                            CHST(iCBC).RxSymEst(:,MCIDX) = xyzh;
                        end
                    case 3
                        % MAP estimate (on tristimulus values)
                        % calibrate for primary powers
                        th(1) = Xh(1)/RGBLED(iCBC).PSDs{1}.rdFlux;
                        th(2) = Xh(2)/RGBLED(iCBC).PSDs{2}.rdFlux;
                        th(3) = Xh(3)/RGBLED(iCBC).PSDs{3}.rdFlux;
                        th = th/sum(th);
                        
                        vYSym = zeros(3,1);
                        for iSym = 1:M
                            vYSym(1) = th(1)-TRIS(1,iSym);
                            vYSym(2) = th(2)-TRIS(2,iSym);
                            vYSym(3) = th(3)-TRIS(3,iSym);
                            dYSym(iSym) = sum((vYSym.*vYSym),1);
                        end
                        if fSAVECHST
                            CHST(iCBC).RxSymEst(:,MCIDX) = th;
                        end
                    otherwise
                        error('Decoder not defined');
                end
                %             SYMIDXh = find(dYSym == min(dYSym),1,'first');
                [~,SYMIDXh] = min(dYSym,[],2);
                if fSAVECHST
                    CHST(iCBC).RxIdx(MCIDX) = SYMIDXh;
                end
                BITSh = dec2binMat(SYMIDXh-1,log2(M));
                
                BITERR = BITERR + biterr2(BITS, BITSh);
                BITCOUNT = BITCOUNT + BITSSYM;
            end
            
            BER(iSNR,iCBC) = BITERR/BITCOUNT;
            if fSAVECHST
                CHST(iCBC).BER = BER(iSNR,iCBC);
                % Determine if state should be saved or not
                FLGST = false;
                if (IDXSNRST(iCBC) <= LENSNRST)
                    if (RNGSNRDB(iSNR,iCBC) > RNGSNRST(IDXSNRST(iCBC)))
                        FLGST = true;
                        IDXSNRST(iCBC) = IDXSNRST(iCBC) + 1;
                    end
                end
                
                if (IDXBERST(iCBC) <= LENBERST)
                    if (BER(iSNR,iCBC) < RNGBERST(IDXBERST(iCBC)))
                        FLGST = true;
                        IDXBERST(iCBC) = IDXBERST(iCBC) + 1;
                    end
                end
                
                IDXCHST(iCBC) = IDXCHST(iCBC) + 1;
                if fSAVECHST
                    FileChnlSt = [ctFileChnlStPRE sprintf('_CBC%d_%d',fCBC,IDXCHST(iCBC)) CHARIDXARCHIVE '.mat'];
                    save(FileChnlSt,'CHST');
                end
            %             clear CHST;
            end
            % Calculate change in BER and SNR
            if iSNR>1
                DSNR = RNGSNRDB(iSNR,iCBC) - RNGSNRDB(iSNR-1,iCBC);
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
            if (RNGSNRDB(iSNR,iCBC) >= RNGSNRMAX) || (BER(iSNR,iCBC) < BERTHMIN)
                LOOPDONE = true;
                LOOPCOUNT = iCBC*RNGSNRLOOP;
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
    end
    % end loop over CSK CBCs
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    %% prep stop
    RNGSNRDB(BER==0) = nan;
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
catch ex
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
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
fprintf('--scrCSK Done--\n');
scrCSKPL;                                                  % Call Show Results script






























































