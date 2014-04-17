close all;
clearvars;
clc;
%% Explore SVD-VLC model
% Generate some random locations in room
% Calculate 'H' matrix
% Apply the same illumination constraint
% Generate pseudorandom sequence
% Calculate illumination
% Calculate capacity

%% FLAGS
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fFILTBLUE = false;
fNOISESHOTONLY = false;

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\3. s-MIMO\5. CapacitySVDKnownCSI\SVD-VLC ex\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrSVDVLC.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/3. s-MIMO/5. CapacitySVDKnownCSI/SVD-VLC ex/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrSVDVLC.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\3. s-MIMO\\5. CapacitySVDKnownCSI\\SVD-VLC ex\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrSVDVLC.m';
    otherwise
        error('Station not defined');
end

ctFileCodeDest = [ctDirRes 'scrSVDVLC.m'];

%%  Set variables
lkIl = 400;     % illumination level at center of room
ctIllTh = 300;
lkPl = 1;       % plane for illumination
txa = 5e-2;     % transmitter side length
txD = 0.5;      % transmitter pitch
txm = 2;        % transmitter m
opf = 5e-3;     % focal length
opD = 5e-3;     % aperture diameter
opFOV = pi/3;   % optics FOV
opfN = opf/opD; % f/#
rxa = 2e-3;     % sensor side length
rxal = 1e-3;    % pixel side length
rxB = 50e6;     % receiver bandwidth
fName = 'Helvetica';
fSize = 16;

%% Configuration
rnSEED = 1;
rng(rnSEED);
NUMFRAMES = 1024;
arrX = 0.2:0.2:3.8;
arrY = 0.2:0.2:3.8;
arrZ = 0:0.2:2.4;
NITER = 1;
% rndX = arrX(random('unid',numel(arrX),NITER,1));
% rndY = arrY(random('unid',numel(arrY),NITER,1));
% rndZ = arrZ(random('unid',numel(arrZ),NITER,1));
rndX = 2;
rndY = 2;
rndZ = 1;
fSize = 16;
fName = 'Helvetica';
CAPLOC = zeros(1,NITER);
%% constants
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX; % wavelengths

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
% Wch1 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.25*Wpsd);   % luminaire SPD
% Wch2 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.25*Wpsd);   % luminaire SPD
% Wch3 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,1*Wpsd);   % luminaire SPD
% Wch4 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.25*Wpsd);   % luminaire SPD

Ambpsd = 5.8e-2*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

if fFILTBLUE        % if blue filtering selected
    Bf = zeros(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
end

%% initialize room
room = cRoom(4,4,4);        % create room and set size
[txX, txY, txZ] = getGrid(2,2,1,txD,txD,2);
% [txX, txY, txZ] = getGrid(room.L,room.W,1,txD,txD,2,'Fill');
txX = txX + room.L/2;
txY = txY + room.W/2;
txZ = txZ + 3;
Ntx = numel(txX);
txLoc = cLocation(txX,txY,txZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% lSP = [1 1 3;3 3 3];
plSP = [2 2 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%
pPratio = ones(Ntx,1);

% for iSP = 1:1:size(lSP,1)
%     idXYlmSP = intersect(find(txX==lSP(iSP,1)),find(txY==lSP(iSP,2)));
%     pPratio(idXYlmSP) = 20;
% end
% txLoc1 = cLocation(txX(1),txY(1),txZ(1));
% txLoc2 = cLocation(txX(2),txY(2),txZ(2));
% txLoc3 = cLocation(txX(3),txY(3),txZ(3));
% txLoc4 = cLocation(txX(4),txY(4),txZ(4));
txOri = cOrientation(pi,0,0);
for iL = 1:1:Ntx
    Wch(iL) = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,pPratio(iL)*Wpsd);   % luminaire SPD
    aLum(iL) = cLuminaire(Wch(iL),cLocation(txX(iL),txY(iL),txZ(iL)));
    aLum(iL).order = txm;
    aLum(iL).orientation = txOri;
    aLum(iL).dimension = cSize(txa,txa,1e-3);
end
room.luminaire = aLum;
% room.setIlluminance(cLocation(txX(idXYlmSP),txY(idXYlmSP),lkPl),lkIl);
room.setIlluminance(cLocation(plSP(1),plSP(2),plSP(3)),lkIl);
%%%%%%
% Pavg = room.luminaire.lmFlux;
pP = zeros(numel(txX),1);     % p-illumination constraint
for iL = 1:1:Ntx
    pP(iL) = room.luminaire(iL).lmFlux;
end

%%%%%%
[plX, plY, plZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
% [plX, plY, plZ] = getGrid(1,1,1,2,2,2);
plX = plX + room.L/2;
plY = plY + room.W/2;
plZ = plZ + 1;
plLoc = cLocation(plX,plY,plZ);
NPloc = numel(plX);
idXYplSP = intersect(find(plX==plSP(1,1)),find(plY==plSP(1,2)));

NPxX = rxa/rxal;
NPxY = rxa/rxal;
NPxD = rxal;
[pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);

if fSAVEALL
    delete([ctDirRes '*.png']);
    delete([ctDirRes '*.mat']);
end

% create buffer to store t-links
tX = zeros(Ntx,NUMFRAMES);

ifig = 1;
for iter=1:1:NITER
    rxX = rndX(iter);
    rxY = rndY(iter);
    rxZ = rndZ(iter);
    rxLoc = cLocation(rxX,rxY,rxZ);
    Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,opf,opfN);
    if fFILTBLUE
        Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
        Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
    end
    Wrx.optics.T = 1; % set lens transmission to 100%
    Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;
    Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;
    Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;
    Wrx.optics.fov = opFOV;
    
    
    % for each receiver group
    disp('Calculating Channel Matrix');
    H = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);
    fsH = reshape(room.getFreeSpaceGain(plLoc,cOrientation(0,0,0),pi/2),numel(plX),sum(room.lmCount(:)));
    
    nL = numel(room.luminaire);
%     txSz = room.lmDimension;
    % currently assuming just 1 lum group at multiple loc
    color = room.luminaire(1).color;
    txOri = room.luminaire(1).orientation;
    Hs = reshape(H,Wrx.rxCount,sum(room.lmCount(:)));
    disp('Calculating Signal and Noise');
    
    % to get channel matrix (A/W), use normalized color (power 1W)
    % % Isig will then be A/W i.e. effective channel matrix H
    [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp, Hch] = ...
        Wrx.getSignal(color,Hs,Ambch,txLoc,cSize(txa,txa,1e-3),txOri);
    % [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = ...
    %     Wrx.getSignal(color.*(color.lmFlux),Hs,Ambch,txLoc,txSz,txOri);
    Hch1 = squeeze(Hch(1,:,:))';
    
    if fNOISESHOTONLY
        In = 0;
    else
        In = Wrx.tia.getNoise(rxB);
    end
    Ish = 2*1.6e-19*(Isig+Iamb)*rxB;
    % Ishc = U'*I2';
    [U, Amt, V] = svd(Hch1);
    A = diag(Amt);
    G = rank(Hch1);
    tP = V'*pP;     % t-illumination constraint
    flPreFix = sprintf('%0.0d_ITER__rxX_%0.1f_rxY_%0.1f_rxZ_%0.1f_Rank_%0.0d__',iter,rxX,rxY,rxZ,G);
    ctFileVars = [ctDirRes flPreFix 'SVDVLCdata.mat'];
    %% Logic
    disp(flPreFix);
    if G>0
        tX(1:G,:) = 2*repmat(tP(1:G),1,NUMFRAMES).*rand(G,NUMFRAMES);
    end
    if G < Ntx
        % generate I1 streams
        tX(G+1:end,:) = repmat(tP(G+1:end),1,NUMFRAMES);
    end
    % Multiplex I1 and I2 streams to transmit over p-space
    pX = V*tX;
    pEX = mean(pX,2);
    pvIll = sum(repmat(pEX',size(fsH,1),1).*fsH,2);
    
    %     for i=1:1:G
    %         CAPLOC(iter) = CAPLOC(iter) + log2(1 + ((Pj(i)*A(i))^2)/In(1));
    %     end
    
    figure;
    surf(plX,plY,reshape(pvIll,size(plX)),'FaceColor','interp');
    minY = min(pvIll(:));
    maxY = max(pvIll(:));
    if maxY == minY,
        maxY = maxY + 1;
    end
    set(gca,'climmode','manual');
    set(gca,'clim',[0 450]);
    cb = colorbar;
    set(cb,'ylim',[0 450], 'ytick',0:50:450,'FontSize',fSize,'FontName',fName);
    axis([0 4 0 4 0 450]);
    set(gca,'FontSize',fSize,'FontName',fName);
    grid on;
%     colorbar;
    view(3);
    rotate3d on;
    text('Position',[plX(idXYplSP) plY(idXYplSP) pvIll(idXYplSP)],'Interpreter','latex','String','{$${\bf{X}}$$}','FontSize',fSize,'FontName','FixedWidth','HorizontalAlignment','center');
    tStr = sprintf('Illuminance\nSetpoint(X): %0.2f lx, Max: %0.2f lx',pvIll(idXYplSP),max(pvIll(:)));
    title(tStr,'FontSize',fSize,'FontName',fName);
    xlabel('Length','FontSize',fSize,'FontName',fName);
    ylabel('Width','FontSize',fSize,'FontName',fName);
    if fSAVEALL
        fname = [ctDirRes 'mIlluminance.fig'];
        saveas(gcf,fname,'mfig');
        fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 3D.png'];
        saveas(gcf,fname,'png');
        view(2);
        fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 2D.png'];
        saveas(gcf,fname,'png');
        ifig = ifig + 1;
%         view(3);
        if fCLOSEFIGS
            close;
        end
    end
    %% FIGURES
   
    % draw spots
    figure('Visible','on');
    h1 = gca;
    Wrx.getImage(cLocation(rxX,rxY,rxZ),txLoc,cSize(txa,txa,1e-3),txOri,h1);
    if fSAVEALL
        fname = [ctDirRes flPreFix num2str(ifig) ' Spots.png'];
        saveas(gcf,fname,'png');
        ifig = ifig + 1;
    end
    if fCLOSEFIGS
        close;
    end
    
    %% SAVE workspace
    disp('Saving Results');
    if fSAVEALL
        save(ctFileVars);
    end
end

    locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
% plot illuminance,irradiance
% figure('Visible','off');
figure;
H(1) = gca;
rotate3d on;
% figure('Visible','off');
figure;
H(2) = gca;
rotate3d on;
room.drawIlluminance(plLoc,ctIllTh,H);
if fSAVEALL
    set(0,'CurrentFigure',get(H(1),'Parent'));
    fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 3D.png'];
    saveas(gcf,fname,'png');
    view(2);
    fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 2D.png'];
    saveas(gcf,fname,'png');
    view(3);
    if fCLOSEFIGS
        close;
    end
    ifig = ifig + 1;
    set(0,'CurrentFigure',get(H(2),'Parent'));
    fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance Coverage.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
    if fCLOSEFIGS
        close;
    end
end

% figure('Visible','off');
figure;
room.drawIrradiance(plLoc);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' Irradiance 3D.png'];
    saveas(gcf,fname,'png');
    view(2);
    fname = [ctDirRes flPreFix num2str(ifig) ' Irradiance 2D.png'];
    saveas(gcf,fname,'png');
    view(3);
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
else
    rotate3d on;
end

% draw setup
% figure('Visible','off');
figure;
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

% % plot LED PSDs
% % figure('Visible','Off');
% figure('Visible','off');
% plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax);
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Normalized PSD');
% % title('Normalized PSDs of LEDs');
% if fSAVEALL
%     fname = [ctDirRes flPreFix num2str(ifig) ' LED PSD.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
% % CIE plot
% % figure('Visible','Off');
% figure('Visible','off');
% obs = cCIE;
% obs.getCoordinates(Wch.npsd);
% if fSAVEALL
%     fname = [ctDirRes flPreFix num2str(ifig) ' CIEplot.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
% % plot filter responses
% % figure('Visible','Off');
% figure('Visible','off');
% plot(Wrx.sensor.filter(1).X,Wrx.sensor.filter(1).Y);
% title('Filter response');
% if fSAVEALL
%     fname = [ctDirRes flPreFix num2str(ifig) ' FilterResponse.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
% % plot receiver responses
% % figure('Visible','Off');
% figure('Visible','off');
% plot(Wrx.sensor.responsivity(1).X,Wrx.sensor.responsivity(1).Y);
% % title('Sensor Responsivity');
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Sensor Responsivity');
% if fSAVEALL
%     fname = [ctDirRes flPreFix num2str(ifig) ' SensorResponsivity.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
% 
% % Draw Cap vs Ill vs Aperture plot
% figure('Visible','off');
% plot(rngDim,C_LID(idLk2,:,1));
% lgd{1} = ['Do = ' num2str(rngOpD(1)*1e3,'%0.2f') 'mm'];
% hold all;
% for idxOpD = 2:size(C_LID,3)
%     opD = rngOpD(idxOpD);
%     plot(rngDim,C_LID(idLk2,:,idxOpD));
%     lgd{idxOpD} = ['Do = ' num2str(rngOpD(idxOpD)*1e3,'%0.2f') 'mm'];
% end
% grid on;
% xlabel('Dimming (%)');
% ylabel('Capacity (b/s/Hz)');
% legend(lgd);
% % axis([min(rngIllSet*100/ctDimRef100) max(rngIllSet*100/ctDimRef100) min(C_LID(:))/1e9 max(C_LID(:))/1e9]);
% axis auto;
% set(gca, 'xdir','reverse');
% if fSAVEALL
%     fname = [ctDirRes 'fig_CapvsDimvsAp.fig'];
%     saveas(gcf,fname,'mfig');
%     fname = [ctDirRes 'CapvsDimvsAp.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
% 
% % Draw Cap vs LinkD vs Aperture plot
% figure('Visible','off');
% plot(rngLkD,C_LID(:,idIll400,1));
% lgd{1} = ['Do = ' num2str(rngOpD(1)*1e3,'%0.2f') 'mm'];
% hold all;
% for idxOpD = 2:size(C_LID,3)
%     opD = rngOpD(idxOpD);
%     plot(rngLkD,C_LID(:,idIll400,idxOpD));
%     lgd{idxOpD} = ['Do = ' num2str(rngOpD(idxOpD)*1e3,'%0.2f') 'mm'];
% end
% grid on;
% xlabel('Link Distance (m)');
% ylabel('Capacity (b/s/Hz)');
% legend(lgd);
% hold off;
% axis auto;
% if fSAVEALL
%     fname = [ctDirRes 'fig_CapvsLinkvsAp.fig'];
%     saveas(gcf,fname,'mfig');
%     fname = [ctDirRes 'CapvsLinkvsAp.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% if fCLOSEFIGS
%     close;
% end
