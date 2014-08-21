% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% FLAGS
fSTATION = 4;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = true;
fCLOSEALL = false;
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

% CIE1931 RGB monochromatic wavelengths
% RMN = 700; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Red led
% GMN = 546; GSC = 1; GSD = 0;               % Mean, SD and scale to generate SPD of Green led
% BMN = 436; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led

RMN = 703; RSC = 1; RSD = 0;              % Mean, SD and scale to generate SPD of Red led
GMN = 564; GSC = 1; GSD = 0;               % Mean, SD and scale to generate SPD of Green led
BMN = 429; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led
% BMN = 509; BSC = 1; BSD = 0;               % Mean, SD and scale to generate SPD of Blue led
cieFile = 'CIE1931_JV_1978_2deg';                 % CIE XYZ CMF curves file
flCIE = [cieFile '.csv'];
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
SPDTYP = SPDTYPE.GAUSSIAN;
cResp = cResponsivity();                    % Responsivity class instance
RResp = cResp.getResponsivity(RMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
GResp = cResp.getResponsivity(GMN);    % Get responsivities vs wavelength for Si-PIN PD (default)
BResp = cResp.getResponsivity(BMN);    % Get responsivities vs wavelength for Si-PIN PD (default)

NTX = 3; NRX = 3;

TOTALBITS = 2e3;                            % Total bit for transmtter to simulate

WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running CSK Simulation...'; % Wait Box title
FIGTITLE = 'Off';
MKTYP = {'o','+','^','s'};
% MKCLR = {'g','y','b','r'};
MKCLR = {[0 1 0],[1 0.5 0],[0 0 1],[1 0 0]};

%% M-CSK CONSTELLATION
M = 4;
%      00,    01,    10,    11     (g w b r)
% CBC 1
x = [0.402; 0.435; 0.169; 0.734];
y = [0.597; 0.290; 0.007; 0.265];
% CBC 2: this is as specified in standard. varies slightly from actual
% calculation. !!                                                           **********IMP**********
% x = [0.402; 0.382; 0.011; 0.734];
% y = [0.597; 0.532; 0.733; 0.265];
Yc = 1;

%% ranges
RNGSNRMIN = 0; RNGSNRMAX = 100; SNROFST = 0;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000]; %DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR
BERTH = 1e-3;   BERTHMIN = 0.75*BERTH;       % BER thresholds;

%% config

% STATION
switch fSTATION
    % Results directory; Spource file; LED table dir;
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\14. CSK\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrCSK.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/14. CSK/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrCSK.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\14. CSK\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSK.m';
    case 4
        ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\14. CSK\\';
        ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSK.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrCSK' CHARIDXARCHIVE '.m']; % Script copy name
ctFileVars = [ctDirRes STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];   % Data file name
if ~exist(ctDirRes,'dir')   % if data directory does NOT exist
    mkdir(ctDirRes);        % create data dir
end

%% logic
try
    %% prep start
    % Wait Bar to show progress
    hWB = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
        'Visible','Off',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
    setappdata(hWB,'canceling',0);
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
    % STATION
    switch fSTATION
        % Results directory; Spource file; LED table dir;
        case 1
            ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\';
            sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
        case 2
            ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/LEDPSD/';
            sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '/' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '/'];
        case 3
            ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
            sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
        case 4
            ctMatDir = 'C:\\Users\\Pankil\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
            sPSDDIR = [ctMatDir cieFile '\' sSPDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
        otherwise
            error('Station not defined');
    end
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
    Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rspd);   % Normalized RED SPD
    Gspd = Gspd/(sum(Gspd)*LAMBDADELTA);
    Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gspd);   % Normalized GREEN SPD
    Bspd = Bspd/(sum(Bspd)*LAMBDADELTA);
    Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bspd);   % Normalized BLUE SPD
    
    %% variables
    if exist(RGBledmat,'file')              % If LED characterization table exists
        load(RGBledmat,'RGB');              % Load RGB led
    else
        RGB = cLEDrgb(RES,Rch,Gch,Bch,flCIE);     % Create RGB led
        RGB.initialize();                   % Initialize led
        if ~exist(sPSDDIR,'dir')
            mkdir(sPSDDIR);                 % Create directory to store characterization data
        end
        save(RGBledmat,'RGB');              % Save LED characterization
    end    
    
    BITERR = 0;
    RNGSNRDB = [];
    
    %% Generate Symbols
    SYMS = zeros(NTX,M);
    PTXAVG = 0;
    for iSym = 1:M
        xc = x(iSym); yc = y(iSym);
        [~,~,~,~,tr,tg,tb] = RGB.getPSD(xc,yc);       % Get transmitter(s) SPD for set x,y
        SYMS(1,iSym) = tr;
        SYMS(2,iSym) = tg;
        SYMS(3,iSym) = tb;
        PTXAVG = PTXAVG + sum(SYMS(:,iSym),1)/4;
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
%         figure('Name','Constellation','NumberTitle',FIGTITLE);
%         RGB.obs.showGamut();
%         hold on;
        title(sprintf('Signal constellation for SNR_{avg}^{tx} = %0.2f(dB)',RNGSNRDB(iSNR)));

        while(BITCOUNT < TOTALBITS - BITSSYM + 1)
            BITS = randi([0 1],[1 BITSSYM]);            % Generate random bits
            SYMIDX = bin2decMat(BITS)+1;
            
            X = SYMS(:,SYMIDX);
            
            % Compute noise
            W = (PTXAVG/SNRrt)*randn(NRX,1);
%             Y = H*X;
            Y = H*X + W;
            
            %% Estimate transmitted vector
            Xh = H\Y;
            XYZh = RGB.Txyz*Xh;
            xyzh = XYZh/sum(XYZh,1);
%             scatter(xyzh(1),xyzh(2),MKTYP{SYMIDX},'MarkerEdgeColor',MKCLR{SYMIDX},'linewidth',2);
           
            % MAP estimate (on RGB data)
%             for iSym = 1:M
%                 vYSym = Y-SYMSRX(:,iSym);
%                 dYSym(iSym) = sum((vYSym.*vYSym),1);
%             end
%             SYMIDXh = find(dYSym == min(dYSym),1,'first');
%             BITSh = dec2binMat(SYMIDXh-1,log2(M));
            
            %% MAP estimate (on xyz data)
            vYSym = zeros(3,1);
            for iSym = 1:M
                vYSym(1) = xyzh(1)-x(iSym);
                vYSym(2) = xyzh(2)-y(iSym);
                dYSym(iSym) = sum((vYSym.*vYSym),1);
            end
            SYMIDXh = find(dYSym == min(dYSym),1,'first');
            
            BITSh = dec2binMat(SYMIDXh-1,log2(M));

            BITERR = BITERR + biterr2(BITS, BITSh);
            BITCOUNT = BITCOUNT + BITSSYM;
        end
%         scatter(x,y,80,'k','x','linewidth',2);
%         grid on;
        
        BER(iSNR) = BITERR/BITCOUNT;
        % Calcculate change in BER and SNR
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
        waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %0.0f min',PROGRESS*100,TREM/60));
        if(getappdata(hWB,'canceling'))
            delete(hWB);
            error('Simulation aborted');
        end
                            
        % check and break if BER for ALL channels are below set thresholds
        if (RNGSNRDB(iSNR) >= RNGSNRMAX) || (BER(iSNR) < BERTHMIN)
            LOOPDONE = true;
            LOOPCOUNT = RNGSNRLOOP;
            PROGRESS = LOOPCOUNT/TOTALLOOPS;
            waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...',PROGRESS*100));
            if(getappdata(hWB,'canceling'))
                delete(hWB);
                error('Simulation aborted');
            end
        end
        iSNR = iSNR + 1;
    end
    %% prep stop
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
    delete(hWB);
catch ex
    delete(hWB);
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






























































