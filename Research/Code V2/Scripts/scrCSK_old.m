% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% config
RNGCBC = 1;                               % CBCs to consider
M = power(2,2);
TOTALBITS = 1e4;                            % Total bit for transmtter to simulate
% DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR

% FLAGS
fSAVEALL = true;
fCLOSEALL = true;
fSAVECHST = true;
fDECODER = 2; % 1.RGB 2.XYZ 3.TRIs
fSHOWPGBAR = isequal(strfind(pwd,'graduate/pbutala'),[]);
fARCHIVE = false;
CHAROVERWRITE = '~';
STRPREFIX = sprintf('M%02d_',M);
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
rng('default');

%% set paths
fs = filesep;
ctFileCodeSrc = [mfilename('fullpath') '.m'];                           % get fullpath of current file
[ctScrDir,~,~] = fileparts(ctFileCodeSrc);                              % get scripts dir
cd(ctScrDir);                                                           % set scripts dir as pwd (reference)
ctDirRes = ['..' fs '..' fs '..' fs '..' fs 'MatlabResults' fs '17. CSK' fs];
ctDirOFDM = ['..' fs '..' fs '..' fs '..' fs 'OFDMcode' fs];
ctDirData = [ctDirRes STRPREFIX 'Data' fs];
% ctDirRes = '../../../../MatlabResults/14. CSK/';
% ctDirData = [ctDirRes STRPREFIX 'Data/'];

ctFileCodeDest = [ctDirData STRPREFIX 'scrCSK' CHARIDXARCHIVE '.m'];    % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];      % Data file name
ctFileChnlStPRE = [ctDirData STRPREFIX 'datChnlStat'];                  % Channel state file name
if ~exist(ctDirData,'dir')                                              % if data directory does NOT exist
    mkdir(ctDirData);                                                   % create data dir
end
ctMatDir = ['..' fs 'Matfiles' fs 'LEDPSD' fs];
% ctMatDir = '../Matfiles/LEDPSD/';

addpath(genpath('..'));
addpath(genpath(ctDirOFDM));

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
    SYMS = zeros(NTX,M,LENCBC);
    LFLX = zeros(M,LENCBC);
    TRIS = zeros(NTX,M,LENCBC);
    
    TSTART = tic;
    % loop over CSK CBCs
    for iCBC = 1:LENCBC
        fCBC = RNGCBC(iCBC);
        
        % CSK object to get CBCs
        csk(iCBC) = cCSK(fCBC,M);                                       % CSK object
        RMN = csk(iCBC).CBC(1).Center; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Band-i (~Red)
        GMN = csk(iCBC).CBC(2).Center; GSC = 1; GSD = 0;      % Mean, SD and scale to generate SPD of Band-j (~Green)
        BMN = csk(iCBC).CBC(3).Center; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Band-k (~Blue)
        
        % Get responsivity for three bands
        RResp(iCBC) = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        GResp(iCBC) = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        BResp(iCBC) = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
        
        [x,y] = csk(iCBC).getSyms();
        
        sPSDDIR = [ctMatDir cieFile fs sSPDTYP fs sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
            RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) fs];
        RGBledmat = [sPSDDIR sprintf('res_%0.5f',RES) '.mat'];                  % LED table mat-file
        
        %% SPDs
        if exist(RGBledmat,'file')              % If LED characterization table exists
            load(RGBledmat,'RGB');              % Load RGB led
        else
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
            RGB = cLEDs(RES,Rch,Gch,Bch,flCIE);
            RGB.initialize(fSHOWPGBAR);                   % Initialize led
            if ~exist(sPSDDIR,'dir')
                mkdir(sPSDDIR);                 % Create directory to store characterization data
            end
            save(RGBledmat,'RGB');              % Save LED characterization
        end
        RGBLED(iCBC) = RGB;
        RGBLED(iCBC).PSDs{1}.name = sprintf('i band (%dnm)',RMN);
        RGBLED(iCBC).PSDs{2}.name = sprintf('j band (%dnm)',GMN);
        RGBLED(iCBC).PSDs{3}.name = sprintf('k band (%dnm)',BMN);
        PPRIM(:,iCBC) = [RGBLED(iCBC).PSDs{1}.rdFlux; RGBLED(iCBC).PSDs{2}.rdFlux; RGBLED(iCBC).PSDs{3}.rdFlux];
        
        BITERR = 0;
        
        %% Generate Symbols
        for iSym = 1:M
            xc = x(iSym); yc = y(iSym);
            [S,Ds,Ts] = RGBLED(iCBC).getPSD(xc,yc);       % Get transmitter(s) SPD for set x,y
            for iTx = 1:NTX
                TRIS(iTx,iSym,iCBC) = Ts(iTx);
                SYMS(iTx,iSym,iCBC) = Ds{iTx}.rdFlux;
                LFLX(iSym,iCBC) = LFLX(iSym,iCBC) + Ds{iTx}.lmFlux;
            end
            TRIS(:,iSym,iCBC) = TXLM*TRIS(:,iSym,iCBC)/LFLX(iSym,iCBC);
            SYMS(:,iSym,iCBC) = TXLM*SYMS(:,iSym,iCBC)/LFLX(iSym,iCBC);
            PTXAVG(iCBC) = PTXAVG(iCBC) + sum(SYMS(:,iSym,iCBC),1)/M;
            PAVGPERLM(iCBC) = PTXAVG(iCBC)/TXLM;
        end
        XAVG = mean(SYMS(:,:,iCBC),2);
        
        %% Compute channel matrix
        H = [RResp(iCBC) 0 0;0 GResp(iCBC) 0; 0 0 BResp(iCBC)];
        
        %% Compute average received signal power
        SIGRXAVG = sqrt(trace(H*XAVG*XAVG'*H'));
        
        %% Compute SYMS at Receiver
        SYMSRX = zeros(NRX,M);
        for iSym = 1:M
            SYMSRX(:,iSym) = H*SYMS(:,iSym,iCBC);
        end
        
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
                        SYMSEST = SYMS(:,:,iCBC);
                    case 2
                        SYMSEST = [x;y;1-x-y];
                    case 3
                        SYMSEST = TRIS(:,:,iCBC);
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
                X = SYMS(:,SYMIDX,iCBC);
                
                % Compute noise
                Wsd = SIGRXAVG/SNRrt;
                W = (Wsd/NRXrt)*randn(NRX,1);
                
                % channel
                Y = H*X + W;
%                 Y = H*X;
                
                % JUST TO CHECK IF -VE VALUES CAUSE SYMS OUTSIDE GAMUT
                Y(Y<0) = 0;
                % NOT USED FOR PAPER
                
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
                        XYZh = RGBLED(iCBC).Txyz*(Xh./PPRIM(:,iCBC));
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
                            vYSym(1) = th(1)-TRIS(1,iSym,iCBC);
                            vYSym(2) = th(2)-TRIS(2,iSym,iCBC);
                            vYSym(3) = th(3)-TRIS(3,iSym,iCBC);
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






























































