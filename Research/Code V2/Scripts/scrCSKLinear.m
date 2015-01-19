% scrCSKLinear
% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% config
RNGCBC = [2 1 7];                               % CBCs to consider
% RNGCBC = 2;
M = power(2,2);
TOTALBITS = 1e4;                            % Total bit for transmtter to simulate
DELTASNR = [0.01 0.05 0.1 2 3 4 5];                % BER ratios to gracefully calculate next SNR
% DELTASNR = [1 2 5 10 10 10 20];                                                   % SNR increment to gracefully calculate next SNR

% FLAGS
fSAVEALL = true;
fCLOSEALL = false;
fSAVECHST = false;
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
[ctScrDir,ctScrFile,ctScrExt] = fileparts(ctFileCodeSrc);                              % get scripts dir
cd(ctScrDir);                                                           % set scripts dir as pwd (reference)
ctDirRes = ['..' fs '..' fs '..' fs '..' fs 'MatlabResults' fs '14. CSKLinearTr' fs];
ctDirData = [ctDirRes STRPREFIX 'Data' fs];
ctDirOFDM = ['..' fs '..' fs '..' fs '..' fs 'OFDMcode' fs];
ctFileCodeDest = [ctDirData STRPREFIX ctScrFile CHARIDXARCHIVE ctScrExt];    % Script copy name
ctFileVars = [ctDirData STRPREFIX 'datCSKLinear' CHARIDXARCHIVE '.mat'];      % Data file name
ctFileChnlStPRE = [ctDirData STRPREFIX 'datChnlStat'];                  % Channel state file name
if ~exist(ctDirData,'dir')                                              % if data directory does NOT exist
    mkdir(ctDirData);                                                   % create data dir
end

addpath(genpath('..'));
addpath(genpath(ctDirOFDM));

%% CONSTANTS
WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = sprintf('%s Running Simulation...',ctScrFile); % Wait Box title

%% ranges
LENCBC = numel(RNGCBC);                     % number of CBCs to consider

RNGSNRMIN = 0; RNGSNRMAX = 20; SNROFST = 0;
RNGSNRMINPL = 0; RNGSNRMAXPL = 18;
RNGSNRLOOP = RNGSNRMAX - RNGSNRMIN + 1;                                         % Number of SNR in each SNR loop
BERRATIOS = [1 5 10 50 100 500 1000];
BERTH = 1e-3;   BERTHMIN = 0.5*BERTH;       % BER thresholds;

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
    
    for iCBC = 1:LENCBC
        fCBC = RNGCBC(iCBC); 
        % CSK object to get CBCs
        csk(iCBC) = cCSK(fCBC,M);                                       % CSK object
        xi = csk(iCBC).CBC(1).x;
        yi = csk(iCBC).CBC(1).y;
        xj = csk(iCBC).CBC(2).x;
        yj = csk(iCBC).CBC(2).y;
        xk = csk(iCBC).CBC(3).x;
        yk = csk(iCBC).CBC(3).y;
        
        % Power to x,y transformation matrix
        Tijk = [xi,xj,xk;yi,yj,yk;1,1,1];
        
        % get [x;y] for all symbols
        [x,y] = csk(iCBC).getSyms();
        Chr = [x;y;ones(1,numel(x))];
        
        % compute Power for all primaries given constraints in IEEE802.15.7
        P(:,:,iCBC) = Tijk\Chr;
        
        % Compute E[sum(Power)]/sqrt(3) for all symbols. This is used to compute
        % noise SD given SNR_{opt}
        PAVG(:,iCBC) = mean(P(:,:,iCBC),2);
        EP(iCBC) = sum(PAVG(:,iCBC),1);
        
        % channel matrix
        H = eye(3);
        % prepare for simulations
        TSTART = tic;
        BITSSYM = log2(M);
        LOOPDONE = false; iSNR = 1;
        vYSym = zeros(2,1);
        
        while ~LOOPDONE
            if iSNR > 1
                dber = BER(iSNR-1,iCBC);   % Smallest BER from all those above the BERTH limit
                % dynamically select next SNR based on how close the current BER is to threshold
                RNGSNRDB(iSNR,iCBC) = RNGSNRDB(iSNR-1,iCBC) + getDeltaSNR(BERTHMIN,dber,BERRATIOS,DELTASNR);
            else
                RNGSNRDB(iSNR,iCBC) = RNGSNRMIN;                 % ITER 1: Smallest SNR from range
            end
            
            SNR = power(10,RNGSNRDB(iSNR,iCBC)/10);
            BITERR = 0; BITCOUNT = 0;
            
            % monte carlo iterations
            MCIDX = 0;
            while(BITCOUNT < TOTALBITS - BITSSYM + 1)
                MCIDX  = MCIDX + 1;                         % increment monte-carlo index
                BITS = randi([0 1],[1 BITSSYM]);            % Generate random bits
                SYMIDX = bin2decMat(BITS)+1;
                
                X = P(:,SYMIDX,iCBC);
                
                % compute noise
                % StdDev Using sin13a formula
                % SD = EP(iCBC)/SNR;
                % StdDev Using SNR = Tr{HXX'H'}/sigma^2 formula
                SD = sqrt(trace(H*PAVG(:,iCBC)*PAVG(:,iCBC)'*H')/SNR);
                
                W = SD*randn(3,1);
                
                % channel (H = eye(3) -> responsivity is considered unity
                % for all three colors)
                Y = H*X+W;
                
                % estimate transmitted vector under constraint H=eye(3) and sum(P)=1;
                Xh = Y/sum(Y,1);
                
                % estimate x,y that was transmitted
                Chrh = Tijk*Xh;
                xh = Chrh(1); yh = Chrh(2);
                for iSym = 1:M
                    vYSym(1) = xh-x(iSym);
                    vYSym(2) = yh-y(iSym);
                    dYSym(iSym) = sum((vYSym.*vYSym),1);
                end
                [~,SYMIDXh] = min(dYSym,[],2);
                
                % estimate transmitted bits
                BITSh = dec2binMat(SYMIDXh-1,log2(M));
                
                % compute bit errors
                BITERR = BITERR + biterr2(BITS, BITSh);
                
                % compute total bits transmitted
                BITCOUNT = BITCOUNT + BITSSYM;
            end % end monte carlo iterations
            
            % compute BER
            BER(iSNR,iCBC) = BITERR/BITCOUNT;
            
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
        end % end loop over SNR values
    end % end loop over CBCs
    
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    %% prep stop
%     RNGSNRDB(BER==0) = nan;
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
catch ex
    if exist('hWB','var') && ishandle(hWB)
        delete(hWB);
    end
    rethrow(ex);
end

fprintf('--scrCSKLinear Done--\n');
scrCSKPLLinear;                                                  % Call Show Results script






























































