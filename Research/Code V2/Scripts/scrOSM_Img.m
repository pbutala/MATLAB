% scrOSM_Img
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
% STATION
switch fSTATION
    case 1
        ctFileCodeSrcNI = '\\ad\eng\users\p\b\pbutala\Documents\MATLAB\Research\Code V2\Scripts\scrOSM_NonImg.m';
    case 2
        ctFileCodeSrcNI = '/home/pbutala/Documents/MATLAB/Research/Code V2/Scripts/scrOSM_NonImg.m';
    case 3
        ctFileCodeSrcNI = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOSM_NonImg.m';
    otherwise
        error('Station not defined');
end
run(ctFileCodeSrcNI);

close all;
clearvars;
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
set(0,'DefaultFigurePaperPosition',[0 0 8 4]);

% FLAGS
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fARCHIVE = true;
fFILTBLUE = false;

rand('seed',0); % seed for Random Number Generator
randn('seed',0); % seed for Random Number Generator

MODSSK = 1;
MODnSM = 2;
MODgSM = 3;
MODeSM = 4;

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = '_varIMG2';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end
CHARIDXARCHIVENI = '_varNI2';

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\8. OSM All 1-M\5. SNR v BER, IMG\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOSM_Img.m';
        ctDirResNI = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\8. OSM All 1-M\5. SNR v BER, NI\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/8. OSM All 1-M/5. SNR v BER, IMG/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOSM_Img.m';
        ctDirResNI = '/home/pbutala/My Documents/MatlabResults/8. OSM All 1-M/5. SNR v BER, NI/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\8. OSM All 1-M\\5. SNR v BER, IMG\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOSM_Img.m';
        ctDirResNI = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\8. OSM All 1-M\\5. SNR v BER, NI\\';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrOSM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'matOSM' CHARIDXARCHIVE '.mat'];      % file to store workspace

% VARIABLES
%%% Setup
lkIl = 400;
lkPl = 1;       % plane for illumination
% txa = 20e-2;     % transmitter side length
txD = 0.5;      % transmitter pitch
txm = 1;        % transmitter m
opf = 5e-3;     % focal length
% opD = 1e-3;     % aperture diameter
opFOV = pi/3;   % optics FOV
rxa = 2e-3;     % sensor side length
rxal = 1e-3;    % pixel side length
opD = sqrt(4*(rxal.*rxal)/pi); % keep D such that aperture area equals pixel area
opfN = opf/opD; % f/#
rxZAT = cOrientation(0,0,0); % receiver orientation
%%% OSM
rngSNRdb = 120:2.5:250;
rngM = power(2,2);
lenM = length(rngM);

rngNt = power(2,2);
lenNt = length(rngNt);

rngMOD = [2 4];
lenMOD = length(rngMOD);

rngBPS = power(2,2:3);
lenBPS = length(rngBPS);

% rngMS = [0.5 1 sqrt(2) 2];
rngMS = 0.5;
lenMS = length(rngMS);
E0 = 1/sqrt(2);
A0 = 1;
D0 = E0+A0/sqrt(2);
Mg0 = D0*rxal*sqrt(2)/txD;
txa = A0*rxal/(Mg0*sqrt(2));

% CONSTANTS
TOTALBITS = 1e5;
BITSTREAM = randi([0 1],[1,TOTALBITS]);
IDXBRK = 0;

% constants
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX; % wavelengths
s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;

Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);        % White PSD
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


iLgd = zeros(lenNt,lenBPS);
for iBPS=1:lenBPS
    STRBPS = sprintf('%d bpS,',rngBPS(iBPS));
    for iNt=1:lenNt
        STRNT = sprintf('Nt %d,',rngNt(iNt));
        fnameNI = [ctDirResNI STRNT STRBPS 'BER vs SNR' CHARIDXARCHIVENI];
        uiopen([fnameNI '.fig'],1);
        figBER(iNt,iBPS) = gcf;
        
        [~,~,~,text_strings] = legend;
        for l=1:numel(text_strings)
            iLgd(iNt,iBPS) = iLgd(iNt,iBPS) + 1;
            clLgd{iLgd(iNt,iBPS),iNt,iBPS} = text_strings{l};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Logic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iNt = 1:lenNt
    Nt = rngNt(iNt);        % get #tx
    STRNT = sprintf('Nt %d,',Nt);
    for iMS = 1:lenMS
        MS = rngMS(iMS);
        STRMS = sprintf('Ms %0.2f,',MS);
        tM = 1/(MS*Mg0);
        opf = 2/(1+tM);
        opfN = opf/opD; % maintain the same aperture area
        STRF = sprintf('MS %0.2f, fcs %0.0fmm,',MS,opf*1e3);
        disp(STRF)
        
        % assume Nt fits on a rect grid
        k = log2(Nt);
        Ntxx = power(2,ceil(k/2));
        Ntxy = power(2,floor(k/2));
        %     k = sqrt(Nt);
        %     Ntxx = ceil(k);
        %     Ntxy = floor(k);
        %%%%%%%%%%%%%%%%%%%%% SETUP ROOM
        % transmitter setup
        [txX, txY, txZ] = getGrid(Ntxx,Ntxy,1,txD,txD,2); % 2x2 array of transmitters
        % [txX, txY, txZ] = getGrid(room.L,room.W,1,txD,txD,2,'Fill'); % 2x2 array of transmitter
        txX = txX + room.L/2;       % add length location offset
        txY = txY + room.W/2;       % add width location offset
        txZ = txZ + 3;              % add height location offset
        txLoc = cLocation(txX,txY,txZ); % tx location object
        txOri = cOrientation(pi,pi/4,0);   % tx orientation
        txSz = cSize(txa,txa,1e-3);     % tx size
        Ntx = numel(txX);           % compute number of transmitters
        for iL = 1:1:Ntx            % need in case transmitters have different output flux
            Wch(iL) = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);   % luminaire SPD
            aLum(iL) = cLuminaire(Wch(iL),cLocation(txX(iL),txY(iL),txZ(iL)));  % create luminaire object
            aLum(iL).order = txm;                       % set lamb order
            aLum(iL).orientation = txOri;               % set orientation
            aLum(iL).dimension = txSz;                  % set size
        end
        room.luminaire = aLum;      % place luminaires in room
        room.setIlluminance(locCntr,lkIl);        % set illuminance at receiver location
        
        % sensor setup
        NPxX = rxa/rxal; NPxY = rxa/rxal;       % number of pixels
        NPxD = rxal; NPx = NPxX*NPxY;                           % pixel pitch
        [pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2); % pixel location in RCS
        
        % receiver setup
        rxX = room.L/2; rxY = room.W/2; rxZ = 1;  % receiver location (center)
        % rxX = room.L; rxY = room.W; rxZ = 1;  % receiver location (corner)
        % rxX = rxX + txD/2; rxY = rxY + txD/2;   % add offset to location for alignment
        
        rxLoc = cLocation(rxX,rxY,rxZ);         %
        Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,opf,opfN);
        if fFILTBLUE
            Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);  % load blue folter
            Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1)); % receiver responsivity
        end
        Wrx.optics.T = 1;                               % set lens transmission to 100%
        Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;   % set pixel lengths
        Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;   % set pixel widths
        Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;   % set pixel heights
        Wrx.optics.fov = opFOV;                         % set FOV
        Wrx.orientation = rxZAT;                        % set orientation
        
        % calculate p-channel matrix
        dStr = sprintf('Calculating H for Nt=%d x Np=%d',Nt,numel(pxX));
        disp(dStr);
        
        tHfs = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);  % free space channel gains
        fsH = reshape(room.getFreeSpaceGain(rxLoc,cOrientation(0,0,0),pi/2),Wrx.rxCount,sum(room.lmCount(:)));
        Hfs = reshape(tHfs,Wrx.rxCount,sum(room.lmCount(:)));   % reshape H to size=[Npx,Ntx]
        
        % calculate channel matrix
        color = room.luminaire(1).color;                    % get color of luminaires
        txOri = room.luminaire(1).orientation;              % get tx orientations
        [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp,tHp] = ...
            Wrx.getSignal(color,Hfs,Ambch,txLoc,cSize(txa,txa,1e-3),txOri);  % get p-channel matrix
        Hp = squeeze(tHp(1,:,:))';                        % reshape p-channel matrix
        if size(Hp,1) == 1
            Hp = Hp';
        end
        %%%%%%%%%%%%%%%%%%%%% END  SETUP ROOM
        Ptxav = room.luminaire(1).rdFlux;
        for iBPS = 1:lenBPS
            BPS = rngBPS(iBPS);
            STRBPS = sprintf('%d bpS,',BPS);
            
            set(0,'CurrentFigure',figBER(iNt,iBPS));
            axBER = gca;
            set(gca,'YScale','log');
            dlinems = get(0,'DefaultLineMarkerSize');
            set(0,'DefaultLineMarkerSize',6);
            hold all;
            for iMOD = 1:lenMOD
                MODSEL = rngMOD(iMOD);
                STRMOD = [char(OSMtype(MODSEL)) ','];
                %             dStr = sprintf('%s: Generating Symbols',char(OSMtype(MODSEL)));
                disp([STRMOD ': Generating Symbols']);
                
                clear tSYMBS SYMBSNT SYMBSMOD optSYM SYMBS;
                iLgd(iNt,iBPS) = iLgd(iNt,iBPS) + 1;
                % create PAM symbols
                % generate all symbols
                switch MODSEL
                    case MODnSM                 % nSM BPS = log2(Nt*M)
                        %                 Nt = numel(txX);
                        M = ceil(power(2,BPS)/Nt);
                        %                     Xs = getPAMsyms(M,1,1);   % get PAM symbols
                        Xs = (1:M)*2/(M+1);
                        clsyms = cell(1,Nt);
                        [clsyms{:}] = deal(Xs(1:end));
                        tsyms = blkdiag(clsyms{:});
                        tSYMBS = full(tsyms);
                        %                         clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SMod: f=%0.0fmm',opf*1e3);
                        %                             clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SMod: \\eta_{s}=%0.2f',ES);
                        clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SM: Img');
                    case MODeSM                 % eSM BPS = Nt*log2(M)
                        %                 Nt = numel(txX);
                        M = power(2,ceil(BPS/Nt));
                        
                        tSYMBS = getExtSMsyms(Nt,M,1);
                        %                         clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SMux: f=%0.0fmm',opf*1e3);
                        %                             clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SMux: \\eta_{s}=%0.2f',ES);
                        clLgd{iLgd(iNt,iBPS),iNt,iBPS} = sprintf('SMP: Img');
                end
                %             dStr = sprintf('%s: Nt=%d, M=%d, b/Hz=%d',char(OSMtype(MODSEL)),Nt,M,BPS);
                dStr = [STRMOD STRNT STRBPS];
                disp(dStr);
                m = log2(M);
                k = log2(Nt);
                % select floor(log2(#syms)) to transfer data
                %         nBPSYM = floor(log2(size(tSYMBS,2)));
                nSYM = power(2,BPS);
                optSYM = getCodeOpt(tSYMBS);
                SYMBS = optSYM(:,1:nSYM);
                
                %             EXv = mean(SYMBS,2);                  % compare for equivalent average powers
                EXtx = mean(mean(SYMBS,2));
                SYMBS = SYMBS*(Ptxav/EXtx);
                EX = mean(mean(SYMBS,2));
                dStr = sprintf('Mean Sym RdFlux = %0.2f W',EX);
                disp(dStr);
                %         bsStr = sprintf('\n%0.0f b/sym',nBPSYM);
                %         clLgd{iLgd} = [clLgd{iLgd} bsStr];
                
                % COMPUTE BER VS SNR
                %             H = eye(Nt);        % TODO: introuce H from imaging receiver setup
                H = Hp;         % just assigning from that computed above and keeping code below the same.
                HSYMB = zeros(numel(pxX),size(SYMBS,2));
                for iSYM = 1:nSYM
                    HSYMB(:,iSYM) = H*SYMBS(:,iSYM);
                end
                
                %             Havg = sqrt(trace(H*H')/rank(H));
                %             EX = Havg*EXtx;
                bit_err = zeros(length(rngSNRdb),1);
                dStr = sprintf('Start Monte Carlo');
                disp(dStr);
                for idb = 1:length(rngSNRdb)
                    vSNRdb = rngSNRdb(idb);
                    vSNR = power(10,vSNRdb/10);
                    vSNRrt = sqrt(vSNR);
                    err = 0;
                    iBits = 1;
                    iDisp = 1;
                    dStr = sprintf('iter txPavg/N0 %0.2f ',vSNRdb);
                    disp(dStr);
                    while(iBits <= TOTALBITS-BPS+1);
                        bitsBin = BITSTREAM(iBits:iBits+BPS-1);
                        bitsDec = bin2decMat(bitsBin);
                        iBits = iBits + BPS;
                        
                        X = SYMBS(:,bitsDec + 1);
                        
                        W = (EX/vSNRrt)*randn(numel(pxX),1);
                        
                        Y = H*X + W;
                        % Y = H*X;
                        
                        DISTYSYM = zeros(1,nSYM);
                        for iSYM = 1:nSYM
                            DISTYSYM(iSYM) = norm(Y-HSYMB(:,iSYM),'fro').^2;
                        end
                        %                     DYSYM = repmat(Y,1,nSYM) - SYMBS;
                        %                     DISTYSYM = sum((DYSYM.^2),1);
                        bitsDecH = find(DISTYSYM == min(DISTYSYM),1);
                        bitsBinH = dec2binMat(bitsDecH-1,BPS);
                        
                        err = err + biterr2(bitsBin,bitsBinH);
                        
                        if iBits > iDisp*TOTALBITS/10
                            dStr = sprintf('iter bits %0.0f',iBits-1);
                            disp(dStr);
                            iDisp = iDisp + 1;
                        end
                    end
                    bit_err(idb) = err/(iBits-1);
                    if bit_err(idb) < 1e-4
                        bit_err(idb+1:end) = NaN;
                        IDXBRK = IDXBRK + 1;
                        if(IDXBRK > 1)
                            IDXBRK = 0;
                            break;
                        end
                    end
                end
                IDXBRK = 0;
                iLC = rem(iMOD+2,lenLC)+1;
                iLS = rem(iMS+2,lenLS)+1;
                iMK = rem(iMOD+2,lenMK)+1;
                plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
                semilogy(axBER,rngSNRdb,bit_err,plStyle);                % plot bit error vs snr
                disp(' ');
                % end
            end % end MOD
            set(0,'CurrentFigure',figBER(iNt,iBPS));
            %         hold off;
            xlabel('{P_{avg}^{tx}}/{N_{0}} (dB)');
            ylabel('BER');
            bsStr = sprintf('Nt=%d; %0.0f b/sym',Nt,BPS);
            title(bsStr);
            grid on;
            axis tight;
            axis([rngSNRdb(1) rngSNRdb(end) 5e-4 1]);
            legend(axBER,clLgd{1:iLgd(iNt,iBPS),iNt,iBPS},'Location','SouthWest');
            set(0,'DefaultLineMarkerSize',dlinems);
            
            if fSAVEALL
                fname = [ctDirRes STRNT STRBPS 'BER vs SNR' CHARIDXARCHIVE];
                %     fname = [ctDirRes ' OSM BER vs SNR']; % Use this to archive
                saveas(figure(figBER(iNt,iBPS)),[fname '.png'],'png');
                saveas(figure(figBER(iNt,iBPS)),[fname '.fig'],'fig');
%                 saveas(figure(figBER(iNt,iBPS)),[fname '.eps'],'eps');
            end
            
            if fCLOSEFIGS
                close;
            end
        end % end BPS
        
        % draw spots
        figSpots = figure;
        axSpots = gca;
        Wrx.getImage(locCntr,txLoc,txSz,txOri,axSpots);
        strT = sprintf('\\mu_{s}=%0.2f',MS);
        title(strT);
        if fSAVEALL
            fname = [ctDirRes STRNT STRMS 'Spots' CHARIDXARCHIVE];
            saveas(figSpots,[fname '.png'],'png');
            saveas(figSpots,[fname '.fig'],'fig');
%             saveas(figSpots,[fname '.eps'],'eps');
        end
        if fCLOSEFIGS
            close;
        end
    end % end Mg
    %         end % end ES
    %     end % end AS
    % draw setup
    figSetup = figure;
    room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
    if fSAVEALL
        fname = [ctDirRes STRNT 'Setup' CHARIDXARCHIVE];
        saveas(figSetup,[fname '.png'],'png');
        saveas(figSetup,[fname '.fig'],'fig');
%         saveas(figSetup,[fname '.eps'],'eps');
    end
    if fCLOSEFIGS
        close;
    else
        rotate3d on;
    end
end % end Nt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. GENERALIZED SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. EXTENDED SPATIAL MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end



% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
% set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);
disp(' ');
disp('--DONE--');






































