% scr Metameric Modulation
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars -except M N;
clc;
%% set paths
fs = filesep;
ctFileCodeSrc = [mfilename('fullpath') '.m'];                           % get fullpath of current file
[ctScrDir,~,~] = fileparts(ctFileCodeSrc);                              % get scripts dir
cd(ctScrDir);                                                           % set scripts dir as pwd (reference)
addpath(genpath('..'));

% config
xp = 1/3; yp = 1/3; % SET POINT TO GENERATE FOR ALL CBCs
M = 2; % possible (M,N) combinations
N = 4;          % (2,4:5), (4,5:7), (8,7)
MM = cMM(M,N);

CBCs = MM.getCBCs();
RNGMMCBC = CBCs(1:10,:);
LENMMCBC = size(RNGMMCBC,1);                     % number of CBCs to consider

if ~isequal(size(RNGMMCBC,2),M)
    error('Number of symbols(%d) is not equal to M(%d)',size(RNGMMCBC,2),M);
end
TOTALBITS = 2e5;                            % Total bit for transmtter to simulate
DELTASNR = [0.01 0.05 0.1 0.5 1 2 3];                % BER ratios to gracefully calculate next SNR
% DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR

% FLAGS
fCLIPY0 = false;
fSAVEALL = true;
fCLOSEALL = true;
fSAVECHST = false;
fSHOWPGBAR = isequal(strfind(pwd,'graduate/pbutala'),[]);
fARCHIVE = true;
CHAROVERWRITE = '~';
STRPREFIX = sprintf('M%d_N%d_a_',M,N);
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
rng('default');

%% set paths
if fCLIPY0
    ctDirRes = ['..' fs '..' fs '..' fs '..' fs 'MatlabResults' fs '30. MM_clip' fs];
else
    ctDirRes = ['..' fs '..' fs '..' fs '..' fs 'MatlabResults' fs '30. MM_noclip' fs];
end
ctDirOFDM = ['..' fs '..' fs '..' fs '..' fs 'OFDMcode' fs];
ctDirData = [ctDirRes STRPREFIX 'Data' fs];
% ctDirRes = '../../../../MatlabResults/14. CSK/';
% ctDirData = [ctDirRes STRPREFIX 'Data/'];

ctFileCodeDest = [ctDirData STRPREFIX 'scrMM' CHARIDXARCHIVE '.m'];    % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datMM' CHARIDXARCHIVE '.mat'];      % Data file name
ctFileChnlStPRE = [ctDirData STRPREFIX 'datChnlStat'];                  % Channel state file name
if ~exist(ctDirData,'dir')                                              % if data directory does NOT exist
    mkdir(ctDirData);                                                   % create data dir
end
ctMatDir = ['..' fs 'Matfiles' fs 'LEDPSD' fs];
% ctMatDir = '../Matfiles/LEDPSD/';
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

WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running CSK Simulation...'; % Wait Box title
FIGTITLE = 'Off';
MKTYP = {'o','+','^','s'};
% MKCLR = {'g','y','b','r'};
MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};

%% ranges
RNGSNRMIN = 0; RNGSNRMAX = 100; SNROFST = 0;
RNGSNRMINPL = 0; RNGSNRMAXPL = 60;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000];
BERTH = 1e-3;   BERTHMIN = 0.5*BERTH;       % BER thresholds;

if fSAVECHST
    RNGSNRST = 10:10:RNGSNRMAX;   LENSNRST = numel(RNGSNRST);
    RNGBERST = [1e-3 5e-4 1e-4 5e-5];     LENBERST = numel(RNGBERST);
    IDXSNRST = ones(LENMMCBC,1);
    IDXBERST = ones(LENMMCBC,1);
end
IDXCHST = zeros(LENMMCBC,1); % count number of states saved per CBC. used in scrCSKPL

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
    cbcs = cCBC();
    PTXAVG = zeros(LENMMCBC,1);
    PAVGPERLM = zeros(LENMMCBC,1);
    LFLX = zeros(M,LENMMCBC);
    mpC2S = -1*ones(3,M,LENMMCBC); % map CBs to SYM
    mpC2Vstr = {};
    %     mpC2V = zeros(cbcs.CBCNT,LENMMCBC); % map CBs to Vector Index
    
    LOOPCOUNT = 1;
    TOTALLOOPS = LENMMCBC*RNGSNRLOOP;
    RNGSNRDB = [];
    
    TSTART = tic;
    % Get number of transmitters and number of receivers
    for iMMCBC = 1:LENMMCBC
        fRXs = zeros(1,cbcs.CBCNT); % 7 color bands
        for iM=1:M % for each CBC
            CBs = cbcs.getColorBands(RNGMMCBC(iMMCBC,iM));
            for idx = 1:3 % for each CB of a CBC
                fRXs(CBs(idx).id+1) = 1; % find Color Band
                mpC2S(idx,iM,iMMCBC) = CBs(idx).id+1;
            end  % end for each CB of a CBC
        end  % end for each CBC
        % #receivers = #receivers = total # color bands for M symbols.
        NTX = sum(fRXs,2);
        NRX = NTX;
        NRXrt = sqrt(NRX);
        
        mpC2V(:,iMMCBC) = sort(unique(mpC2S(:,:,iMMCBC)));
        mpC2Vstr{iMMCBC} = num2str(mpC2V(:,iMMCBC)'-1,'%d');
        
        SYMS = zeros(NTX,M);
        TRIS = zeros(NTX,M);
        Resp = zeros(NRX,M);
        
        % for each CBC, create RGB LED
        for iM=1:M % for each CBC
            CBs = cbcs.getColorBands(RNGMMCBC(iMMCBC,iM));
            RMN = CBs(1).Center; RSC = 1; RSD = 0;             % Mean, SD and scale to generate SPD of Band-i (~Red)
            GMN = CBs(2).Center; GSC = 1; GSD = 0;             % Mean, SD and scale to generate SPD of Band-j (~Green)
            BMN = CBs(3).Center; BSC = 1; BSD = 0;             % Mean, SD and scale to generate SPD of Band-k (~Blue)
            %             Resp(mpC2S(1,iM,iMMCBC),iM) = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
            %             Resp(mpC2S(2,iM,iMMCBC),iM) = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
            %             Resp(mpC2S(3,iM,iMMCBC),iM) = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
            
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
%             RGBLED(iMMCBC,iM) = RGB;
%             RGBLED(iMMCBC,iM).PSDs{1}.name = sprintf('CB%d (%dnm)',CBs(1).id,RMN);
%             RGBLED(iMMCBC,iM).PSDs{2}.name = sprintf('CB%d (%dnm)',CBs(2).id,GMN);
%             RGBLED(iMMCBC,iM).PSDs{3}.name = sprintf('CB%d (%dnm)',CBs(3).id,BMN);
%             
%             [S,Ds,Ts] = RGBLED(iMMCBC,iM).getPSD(xp,yp);       % Get transmitter(s) SPD for set x,y
              RGBLED = RGB;
              [S,Ds,Ts] = RGBLED.getPSD(xp,yp);       % Get transmitter(s) SPD for set x,y
            
            for iTx = 1:3
                TRIS(find(mpC2V(:,iMMCBC)==mpC2S(iTx,iM,iMMCBC)),iM) = Ts(iTx);
                SYMS(find(mpC2V(:,iMMCBC)==mpC2S(iTx,iM,iMMCBC)),iM) = Ds{iTx}.rdFlux;
                LFLX(iM,iMMCBC) = LFLX(iM,iMMCBC) + Ds{iTx}.lmFlux;
            end
            TRIS(:,iM) = TXLM*TRIS(:,iM)/LFLX(iM,iMMCBC);
            SYMS(:,iM) = TXLM*SYMS(:,iM)/LFLX(iM,iMMCBC);
            PTXAVG(iMMCBC) = PTXAVG(iMMCBC) + sum(SYMS(:,iM),1)/M;
        end % end for each CBC (iM)
        PAVGPERLM(iMMCBC) = PTXAVG(iMMCBC)/TXLM;
        
        BITERR = 0;
        XAVG = mean(SYMS,2);
        
        %% Compute channel matrix
        H = zeros(NRX,NTX);
        % assuming NRX = NTX
        for iTx = 1:NTX
            H(iTx,iTx) = cResp.getResponsivity(cbcs.CBs(mpC2V(iTx,iMMCBC)).Center);
        end
        
        %% Compute average received signal power
        SIGRXAVG = sqrt(trace(H*XAVG*XAVG'*H'));
        
        %% Compute SYMS at Receiver
        SYMSRX = zeros(NRX,M);
        for iSym = 1:M
            SYMSRX(:,iSym) = H*SYMS(:,iSym);
        end
        
        BITSSYM = log2(M);
        LOOPDONE = false; iSNR = 1;
        
        while ~LOOPDONE                                                 % LOOP START SNR (dynamically select next SNR)
            if iSNR > 1
                dber = BER(iSNR-1,iMMCBC);   % Smallest BER from all those above the BERTH limit
                % dynamically select next SNR based on how close the current BER is to threshold
                RNGSNRDB(iSNR,iMMCBC) = RNGSNRDB(iSNR-1,iMMCBC) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
            else
                RNGSNRDB(iSNR,iMMCBC) = RNGSNRMIN;                 % ITER 1: Smallest SNR from range
            end
            
            SNR = power(10,RNGSNRDB(iSNR,iMMCBC)/10);
            SNRrt = sqrt(SNR);
            BITERR = 0; BITCOUNT = 0;
            
            if fSAVECHST
                ARLEN = floor(TOTALBITS/BITSSYM);
                SZRXSYMS = [NRX ARLEN];
                switch fDECODER
                    case 1
                        SYMSEST = SYMS(:,:,iMMCBC);
                    case 2
                        SYMSEST = [x;y;1-x-y];
                    case 3
                        SYMSEST = TRIS(:,:,iMMCBC);
                end
                CHST(iMMCBC) = cChnlState(ARLEN,SZRXSYMS,SYMSEST,H);
                CHST(iMMCBC).SNRdB = RNGSNRDB(iSNR,iMMCBC);
            end
            MCIDX = 0;
            
            while(BITCOUNT < TOTALBITS - BITSSYM + 1)
                MCIDX  = MCIDX + 1;                         % increment monte-carlo index
                BITS = randi([0 1],[1 BITSSYM]);            % Generate random bits
                SYMIDX = bin2decMat(BITS)+1;
                if fSAVECHST
                    CHST(iMMCBC).TxIdx(MCIDX) = SYMIDX;
                end
                X = SYMS(:,SYMIDX);
                
                % Compute noise
                Wsd = SIGRXAVG/SNRrt;
                W = (Wsd/NRXrt)*randn(NRX,1);
                
                % channel
                Y = H*X + W;
                % Y = H*X;
                % clip received symbols
                if fCLIPY0
                    Y(Y<0) = 0;
                end
                
                %% Estimate transmitted vector
                for iSym = 1:M
                    vYSym = Y-SYMSRX(:,iSym);
                    dYSym(iSym) = sum((vYSym.*vYSym),1);
                end
                if fSAVECHST
                    CHST(iMMCBC).RxSymEst(:,MCIDX) = Xh;
                end
                
                [~,SYMIDXh] = min(dYSym,[],2);
                if fSAVECHST
                    CHST(iMMCBC).RxIdx(MCIDX) = SYMIDXh;
                end
                BITSh = dec2binMat(SYMIDXh-1,log2(M));
                
                BITERR = BITERR + biterr2(BITS, BITSh);
                BITCOUNT = BITCOUNT + BITSSYM;
            end % end SNR iteration
            
            BER(iSNR,iMMCBC) = BITERR/BITCOUNT;
            if fSAVECHST
                CHST(iMMCBC).BER = BER(iSNR,iMMCBC);
                % Determine if state should be saved or not
                FLGST = false;
                if (IDXSNRST(iMMCBC) <= LENSNRST)
                    if (RNGSNRDB(iSNR,iMMCBC) > RNGSNRST(IDXSNRST(iMMCBC)))
                        FLGST = true;
                        IDXSNRST(iMMCBC) = IDXSNRST(iMMCBC) + 1;
                    end
                end
                
                if (IDXBERST(iMMCBC) <= LENBERST)
                    if (BER(iSNR,iMMCBC) < RNGBERST(IDXBERST(iMMCBC)))
                        FLGST = true;
                        IDXBERST(iMMCBC) = IDXBERST(iMMCBC) + 1;
                    end
                end
                
                if fSAVECHST
                    if FLGST
                        IDXCHST(iMMCBC) = IDXCHST(iMMCBC) + 1;
                        FileChnlSt = [ctFileChnlStPRE sprintf('_CBC%d_%d',fCBC,IDXCHST(iMMCBC)) CHARIDXARCHIVE '.mat'];
                        save(FileChnlSt,'CHST');
                    end
                end
                clear CHST;
            end
            % Calculate change in BER and SNR
            if iSNR>1
                DSNR = RNGSNRDB(iSNR,iMMCBC) - RNGSNRDB(iSNR-1,iMMCBC);
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
            if (RNGSNRDB(iSNR,iMMCBC) >= RNGSNRMAX) || (BER(iSNR,iMMCBC) < BERTHMIN)
                LOOPDONE = true;
                LOOPCOUNT = iMMCBC*RNGSNRLOOP;
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
        end % end while ~LOOPDONE
    end % end for all MM configs (iMMCBC)
    
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
fprintf('--scrMM Done--\n');
% scrMMPL;                                                  % Call Show Results script










































