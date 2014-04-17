% script OSMimgSVD
close all;
clearvars;
clc;

% Generate different different locations in room.
% for each calculate channel matrix
% calculate SVD
% calculate rank (equal to number of independent transmitters)
% assign PAM/txIndex
% calculate BER vs SNR

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fFILTBLUE = false;
fNOISESHOTONLY = false;

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\5. OSM-SVD3\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOSMimoSVD.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/5. OSM-SVD3/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOSMimoSVD.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\5. OSM-SVD3\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOSMimoSVD.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrOSMimgSVD.m'];

% VARIABLES
% rnglkIl = 50:50:400;     % illumination level at center of room
rnglkIl = 50:50:200;     % illumination level at center of room
rngm = 2:2:8;  % range of #bits/symbol
% rngopD = [1e-3 2e-3 5e-3];
lkPl = 1;       % plane for illumination
txa = 20e-3;     % transmitter side length
txD = 0.5;      % transmitter pitch
txm = 2;        % transmitter m
opf = 5e-3;     % focal length
opD = 1e-3;     % aperture diameter
opFOV = pi/3;   % optics FOV
opfN = opf/opD; % f/#
rxa = 2e-3;     % sensor side length
rxal = 1e-3;    % pixel side length
rxB = 40e6;     % rx bandwidth

% constants
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX; % wavelengths
s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;

Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);        % White PSD
Ambpsd = 5.8e-2*ones(size(lambdas));        % Ambient PSD
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);   % Ambient Channel

% etc
ifig = 1;

if fFILTBLUE        % if blue filtering selected
    Bf = zeros(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
end

% initialize room
room = cRoom(4,4,4);        % create room and set size (4m x 4m x 4m)
locCntr = cLocation(room.L/2,room.W/2,1);

% transmitter setup
% [txX, txY, txZ] = getGrid(2,2,1,txD,txD,2); % 2x2 array of transmitters
[txX, txY, txZ] = getGrid(room.L,room.W,1,txD,txD,2,'Fill'); % 2x2 array of transmitter
txX = txX + room.L/2;       % add length location offset
txY = txY + room.W/2;       % add width location offset
txZ = txZ + 3;              % add height location offset
txLoc = cLocation(txX,txY,txZ); % tx location object
txOri = cOrientation(pi,0,0);   % tx orientation
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

% sensor setup
NPxX = rxa/rxal; NPxY = rxa/rxal;       % number of pixels
NPxD = rxal; NPx = NPxX*NPxY;                           % pixel pitch
[pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2); % pixel location in RCS

% receiver setup
rxX = room.L/2; rxY = room.W/2; rxZ = 1;  % receiver location (center)
% rxX = room.L; rxY = room.W; rxZ = 1;  % receiver location (corner)
rxX = rxX + txD/2; rxY = rxY + txD/2;   % add offset to location for alignment

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

% calculate p-channel matrix
disp('Calculating Channel Matrix');
tHfs = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);  % free space channel gains
fsH = reshape(room.getFreeSpaceGain(rxLoc,cOrientation(0,0,0),pi/2),Wrx.rxCount,sum(room.lmCount(:)));
Hfs = reshape(tHfs,Wrx.rxCount,sum(room.lmCount(:)));   % reshape H to size=[Npx,Ntx]

% calculate channel matrix
color = room.luminaire(1).color;                    % get color of luminaires
txOri = room.luminaire(1).orientation;              % get tx orientations
[Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp,tHp] = ...
    Wrx.getSignal(color,Hfs,Ambch,txLoc,cSize(txa,txa,1e-3),txOri);  % get p-channel matrix
Hp = squeeze(tHp(1,:,:))';                        % reshape p-channel matrix

% calculate noise
if fNOISESHOTONLY
    In = 0;
else
    In = Wrx.tia.getNoise(rxB);                 % calc TIA noise
end
Ish = 2*1.6e-19*(Isig+Iamb)*rxB;                % calc shot noise
I2 = In(1) + Ish;                               % calc total noise (AWGN)

% SVD
[U, Ht, V] = svd(Hp);       % singular value decomposition
G = rank(Ht);               % rank of matrix
Htt = Ht';

% files
flPreFix = sprintf('OSM-SVD_rxX_%0.1f_rxY_%0.1f_rxZ_%0.1f_Rank_%d_',rxX,rxY,rxZ,G);% file prefix to store workspace
disp(flPreFix);                                     % display prefix
ctFileVars = [ctDirRes flPreFix 'OSM-SVD.mat'];      % file to store workspace


% MONTE-CARLO RUNS
bit_err = zeros(numel(rngm),numel(rnglkIl));
Ill = zeros(numel(rngm),numel(rnglkIl));
figSNR = figure;
% set(gca,'YScale','log');
xlabel('Illumination (lx)');
ylabel('Bit Error Ratio');
tStr = sprintf('Imaging OSM');
title(tStr);
grid on;
hold all;
idm = 1;
for m = rngm
    display(sprintf('m = %d', m));
    loop =1;
    for lkIl = rnglkIl
        display(sprintf('Ill target = %0.2f lux', lkIl));
        room.setIlluminance(locCntr,lkIl);        % set illuminance at receiver location
        % calculate p-illumination constraint
        pP = zeros(numel(txX),1);     % p-illumination constraint
        for iL = 1:1:Ntx
            pP(iL) = room.luminaire(iL).lmFlux;
        end
        
        tP = V'*pP;                 % t-illumination constraint
        
        % OSM setup
        Nt = 2^(floor(log2(G)));                   % # of Txs = rank of channel matrix
%         m = 8;                    % bits per symbols
        M = 2^m;                % # of symbols
        mi = 0.9;              % modulation index (provides some power to transmit 0 and successfully detect LED)
        mic = 1-mi;
%         SymsPAM = Nt*(tP(1:Nt)*(mic:mi*2/(M-1):2-mic))'; % Generate PAM symbols
        SymsPAM = Nt*(tP(1:Nt)*(getPAMsyms(M,mi*2/(M-1),1)' + mic))';
        
        Pmt = zeros(Ntx,1);
        err = 0;                % Detected Errors
        iter = 1e3;            % # of Iterations
        for ii = 1:iter
            bit_stream = randi([0,1], [1,log2(M*Nt)]);                  % generate bit stream e.g. 4 bits per symbol for 4 Txs
            sym_index = bin2dec(int2str(bit_stream(1:log2(M))));                % Symbol index 1:log2(M) (begin)
            led_index = bin2dec(int2str(bit_stream(log2(M)+1:log2(Nt*M))));   % transmitter index log2(M)+1:log2(Nt*M) (end)
            tx = zeros(Ntx,1);
            % tx(led_index) = modulate(alphabit,sym_index);             % assign symbol to channel
            tx(led_index+1) = SymsPAM(sym_index+1,led_index+1);                            % assign symbol to channel
            tx(Nt+1:end) = tP(Nt+1:end);
            
            %%%%%%%Add AWGN Channel%%%%%%%%%%%%%%%%%%%%%%%%%%
            % y = awgn((Htt*tx),jj,'measured');         % Y' = HX' + W'
            pW = sqrt(I2(1))*randn(NPx,1);          % get random noise vector
            
            pX = V*tx;
            Pmt = Pmt + pX;
            
            pY = Hp*pX + pW;
            y = U'*pY;
            % y = Ht*tx + pW;                     % add noise to signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%Optimal Spatial Modulation Detector%%%%%%%%%%%%%%
            for j=1:Nt                                                  % antenna (j)
                for q=1:M                                               % symbol (q)
                    % g_jq = Htt(:,i)*modulate(alphabit,j-1);
                    g_jq = Ht(:,j)*SymsPAM(q,j);                          % g_jq = hj*xq
                    min_2(q,j) = norm(g_jq, 'fro')^2 - 2*real(y'*g_jq); % This vector contains the errors to all antennas.
                end
            end
            [tmp,idx] =  min( min_2(:));            % do the argmin{} over all j,q
            [idq,idj] = find(min_2==tmp);         % jointly detect tx index and symbol index
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%Estimated bits from LED index%%%%%%%%%%%%%%%%%%%
            tmp1 = dec2bin(idj-1,log2(Nt));        % generate at least log2(Nt) long binary number.
            % tmp1 has binary index of the estimated LED on. receover bits
            % encoded in LED index.
            for uu = 1:log2(Nt)
                led_bits(uu)= bin2dec(tmp1(uu)); % convert binary string to array of decimal 0/1 values.
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%Estimated bits from Symbol%%%%%%%%%%%%%%%%%%%%%%
            tmp2 = dec2bin(idq-1,log2(M));         % generate at least log2(M) long binary number.
            % tmp2 has binary index of the estimated data. recover bits
            % encoded in data.
            for xx = 1:log2(M)
                sym_bits(xx)= bin2dec(tmp2(xx));    % convert binary string to array of decimal 0/1 values.
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            rx_bits = [sym_bits led_bits];          % concatenate SYM bits and LED bits to recover frame.
            
            err = err + biterr(rx_bits,bit_stream);
        end
        
        pEX(:,loop) = mean(Pmt,2)/iter;
        Ill(idm,loop) = sum(repmat(pEX(:,loop)',size(fsH,1),1).*fsH,2);
        display(sprintf('Ill achieved = %0.02f lux', Ill(idm,loop)));
        bit_err(idm,loop) = err/((log2(M*Nt))*iter);
        display(sprintf('BER = %0.4f', bit_err(idm,loop)));
        loop = loop+1;
        display(sprintf('\n'));
    end
    % - BIT ERROR VS SNR
    semilogy(rnglkIl,bit_err(idm,:));                % plot bit error vs snr
    lgd{idm} = [num2str(M,'%d') '-ary PAM'];
    idm = idm + 1;
end

% delete old files
if fSAVEALL
    delete([ctDirRes '*.png']);
    delete([ctDirRes '*.mat']);
    delete([ctDirRes '*.fig']);
end

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end

figure(figSNR);
legend(lgd);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' BER vs SNR.fig'];
    saveas(gcf,fname,'fig');
    fname = [ctDirRes flPreFix num2str(ifig) ' BER vs SNR.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% PLOT RESULTS

% - draw setup
figSetup = figure;
room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' Setup.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
else
    rotate3d on;
end

% draw spots
figSpots = figure;
h1 = gca;
Wrx.getImage(rxLoc,txLoc,txSz,txOri,h1);
%             set(gca,'FontName',fName,'FontSize',fSize);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' Spots.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% - plot LED PSDs
figPSD = figure;
plot(Wch(1).npsd.X,Wch(1).npsd.Y/Wch(1).npsd.Ymax);
axis tight;
grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized PSD');
% title('Normalized PSDs of LEDs');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' LED PSD.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% - CIE plot
figCIE = figure;
obs = cCIE;
obs.getCoordinates(Wch(1).npsd);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' CIEplot.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% - plot filter responses
figFilt = figure;
plot(Wrx.sensor.filter(1).X,Wrx.sensor.filter(1).Y);
axis tight;
title('Filter response');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' FilterResponse.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% - plot receiver responses
figResp = figure;
plot(Wrx.sensor.responsivity(1).X,Wrx.sensor.responsivity(1).Y);
axis tight;
grid on;
xlabel('Wavelength (nm)');
ylabel('Sensor Responsivity');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' SensorResponsivity.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

figure(figSNR);
figure(figSpots);

% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);




















