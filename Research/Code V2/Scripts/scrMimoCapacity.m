% On capacity of MIMO VLC links with imaging receiver
close all;
clearvars;
clc;

daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');
%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;
fFILTBLUE = false;
fNOISESHOTONLY = false;

%% VARIABLE PARAMS
% rngIllSet = 50:50:800; % varILLSet
rngIllSet = 400; % varILLSet
rngDim = ((max(rngIllSet)-rngIllSet)*100)./max(rngIllSet);
idIll400 = find(rngIllSet == 400);
ctDimRef100 = 400;
% rngIllSet = 400; % varILLSet
% rngTxD = [0.5 0.8 1.0 1.5 2];
% rngLkD = 1:0.2:3;
% rngLkD = 0.2:0.2:3;
% rngLkD = [0.9 1.2 1.9 2.2];
rngLkD = 2;
idLk2 = find(rngLkD == 2);
% rngTxa = [1e-3 1e-2 10e-2 20e-2 50e-2];
% rngOpf = [1e-3 5e-3 1e-2 5e-2 10e-2];
% rngOpD = [1e-3 5e-3 1e-2];
% rngOpD = [1e-3 2e-3 (5e-3*2/sqrt(pi))];
rngOpD = 5e-3;
% rngOpD = [1e-3 5e-3];
% rngRxa = [5e-3 1e-2 5e-2 10e-2 20e-2];
% rngRxal = [0.5e-3 1e-3 5e-3 1e-2 5e-2 10e-2];% depends on rngRxa only use values smaller than 'varRxa'
rotD = 15;
rngRxZeta = 0:rotD*pi/180:pi/2-rotD*pi/180;
% rngRxAlpha = -pi:15*pi/180:pi;
rngRxAlpha = 0:rotD*pi/180:pi-rotD*pi/180;

%%  Set variables
% varIllSet = 400; % varILLSet
% varLkD = 2;
varTxa = 5e-2;
varOpf = 5e-3;
% varOpD = 5e-3;
varOpFOV = pi/3;
varRxa = 5e-3;
varRxal = 1e-3;% depends on rngRxa only use values smaller than 'varRxa'
varTxD = 0.5; %D between TX
% varTxD = sqrt(2)*varRxal*(varLkD-varOpf)/varOpf;
% varRxZeta = pi/24;
% varRxAlpha = 0;
varRxTau = 0;

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\3. s-MIMO\5. CapacitySVDKnownCSI\SVD-VLC full angle\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrMimoCapacity.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/3. s-MIMO/5. CapacitySVDKnownCSI/SVD-VLC full angle/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrMimoCapacity.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\3. s-MIMO\\5. CapacitySVDKnownCSI\\SVD-VLC full angle\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrMimoCapacity.m';
    otherwise
        error('Station not defined');
end


ctFileCodeDest = [ctDirRes 'scrMimoCapacity.m'];
ctIllTh = 300;
ctrxB = 50e6;
%% constants
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
% room.setIlluminance(.) will set the output power of the transmitters

Ambpsd = 5.8e-2*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);
if fFILTBLUE
    %     Bf = getSOG(470,10,1,lambdas);
    Bf = zeros(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
    %     Bf = 1*Bf/max(Bf); % Ideal Blue Filter
end


%% initialize room
room = cRoom(4,4,4);
% [txX, txY, txZ] = getGrid(2,2,1,varTxD,varTxD,2);
[txX, txY, txZ] = getGrid(room.L,room.W,1,varTxD,varTxD,2,'Fill');
% [txX, txY, txZ] = getGrid(1,1,1,2,2,2);
txX = txX + room.L/2;
txY = txY + room.W/2;
txZ = txZ + 3;
txLoc = cLocation(txX,txY,txZ);
txOri = cOrientation(pi,0,0);
room.luminaire = cLuminaire(Wch,txLoc);
room.luminaire.orientation = txOri;
room.luminaire.dimension = cSize(varTxa,varTxa,1e-3);

[plX, plY, plZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
plX = plX + room.L/2;
plY = plY + room.W/2;
plZ = plZ + 1;
plLoc = cLocation(plX,plY,plZ);

NPxX = varRxa/varRxal;
NPxY = varRxa/varRxal;
NPxD = varRxal;
% NPxD = 5.91e-2/NPxX;
[pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);

C_ZA = zeros(numel(rngRxZeta),numel(rngRxAlpha));
ifig = 1;
if fSAVEALL
    delete([ctDirRes '*.png']);
    delete([ctDirRes '*.mat']);
end
for idOpD = 1:numel(rngOpD)
    varOpD = rngOpD(idOpD);
    fN = varOpf/varOpD;
    Wrx.optics.fN = fN;
    for idIllSet = 1:numel(rngIllSet)
        varIllSet = rngIllSet(idIllSet);
        room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
        for idLkD = 1:numel(rngLkD)
            varLkD = rngLkD(idLkD);
            for idZeta = 1:numel(rngRxZeta)
                varRxZeta = rngRxZeta(idZeta);
                for idAlpha = 1:numel(rngRxAlpha)
                    varRxAlpha = rngRxAlpha(idAlpha);
                    
%                     [rxX,rxY,rxZ] = getGrid(1,1,1,2,2,2,'Fill');
                    % [rxX,rxY,rxZ] = getGrid(room.L,room.W,1,1,1,2,'Fill');
                    [rxX,rxY,rxZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
                    rxX = rxX + room.L/2;
                    rxY = rxY + room.W/2;
                    rxZ = rxZ + 3 - varLkD;
                    rxLoc = cLocation(rxX,rxY,rxZ);
                    rxOri = cOrientation(varRxZeta,varRxAlpha,varRxTau);
                    Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,varOpf,fN);
                    if fFILTBLUE
                        Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
                        Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
                    end
                    Wrx.orientation = rxOri;
                    Wrx.optics.T = 1; % set lens transmission to 100%
                    Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;
                    Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;
                    Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;
                    Wrx.optics.fov = varOpFOV;
                    Wrx.sensor.color = [1 1 1];
                    
                    flPreFix = sprintf('%0.0flx_TxD_%0.1fcm_LkD_%0.1fm_Txa_%0.0fmm_Opf_%0.0fmm_OpD_%0.0fmm_Rxa_%0.0fmm_Z_%0.1f_A_%0.1f_T_%0.1f_',...
                        varIllSet,varTxD*1e2,varLkD,varTxa*1e3,varOpf*1e3,varOpD*1e3,varRxa*1e3,varRxZeta*180/pi,varRxAlpha*180/pi,varRxTau*180/pi);
                    ctFileVars = [ctDirRes flPreFix 'MimoCapacity.mat'];
                    %% Logic
                    disp(flPreFix);
                    
                    % for each receiver group
                    disp('Calculating Channel Matrix');
                    H = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);
                    
                    nL = numel(room.luminaire);
                    txSz = room.lmDimension;
                    % currently assuming just 1 lum group at multiple loc
                    color = room.luminaire(1).color;
                    txOri = room.luminaire(1).orientation;
                    Hs = reshape(H,Wrx.rxCount,sum(room.lmCount(:)));
                    disp('Calculating Signal and Noise');
                    
                    % to get channel matrix (A/W), use normalized color (power 1W)
                    % % Isig will then be A/W i.e. effective channel matrix H
                    [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp, Hch] = ...
                        Wrx.getSignal(color,Hs,Ambch,txLoc,txSz,txOri);
                    % [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = ...
                    %     Wrx.getSignal(color.*(color.lmFlux),Hs,Ambch,txLoc,txSz,txOri);
                    Hch1 = squeeze(Hch(1,:,:))';
                    
                    [U, Amt, V] = svd(Hch1);
                    A = diag(Amt);
                    
                    if fNOISESHOTONLY
                        In = 0;
                    else
                        In = Wrx.tia.getNoise(ctrxB);
                    end
                    Ish = 2*1.6e-19*(Isig+Iamb)*ctrxB;
                    I2 = Ish + In(1);
                    % Ishc = U'*I2';
                    
                    Pavg = room.luminaire.rdFlux;
                    
                    % IGNORING SHOT NOISE BECAUSE AMP NOISE DOMINATES
                    %             [L,tP,Pj] = getWaterfill(Pavg, In(1), A);
                    Pj = V'*(Pavg*ones(numel(txX),1));
                    
                    for i=1:1:numel(A)
                        C_ZA(idZeta,idAlpha) = C_ZA(idZeta,idAlpha) + log2(1 + ((Pj(i)*A(i))^2)/In(1));
                    end
                    
                    locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
                    
                    %% SAVE workspace
                    disp('Saving Results');
                    if fSAVEALL
                        save(ctFileVars);
                    end
                    
                    %% FIGURES
                    %                 ifig = 1;
                    
                    % draw spots
                    figure;
                    h1 = gca;
                    Wrx.getImage(locCntr,txLoc,txSz,txOri,h1);
                    %             set(gca,'FontName',fName,'FontSize',fSize);
                    if fSAVEALL
                        fname = [ctDirRes num2str(ifig) '_' flPreFix 'Spots.png'];
                        saveas(gcf,fname,'png');
                        ifig = ifig + 1;
                    end
                    if fCLOSEFIGS
                        close;
                    end
                end
            end
        end
        % plot illuminance,irradiance
        H(1) = gca;
        rotate3d on;
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
    end
end
% draw setup
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

% plot LED PSDs
figure;
plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax);
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
% CIE plot
figure;
obs = cCIE;
obs.getCoordinates(Wch.npsd);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' CIEplot.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% plot filter responses
figure;
plot(Wrx.sensor.filter(1).X,Wrx.sensor.filter(1).Y);
title('Filter response');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' FilterResponse.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% plot receiver responses
figure;
plot(Wrx.sensor.responsivity(1).X,Wrx.sensor.responsivity(1).Y);
% title('Sensor Responsivity');
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

% Draw Cap vs zeta 
figure;
plot(180*rngRxZeta/pi,C_ZA(:,1)*180/pi);

grid on;
xlabel('Zenith Angle (deg)');
ylabel('Capacity (b/s/Hz)');
axis auto;
if fSAVEALL
    fname = [ctDirRes 'fig_CapvsZenith.fig'];
    saveas(gcf,fname,'mfig');
    fname = [ctDirRes 'CapvsZenith.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% Draw Cap vs Azimuth vs Zenith
figure;
plot(180*rngRxAlpha/pi,C_ZA(1,:));
lgd{1} = ['\zeta = ' num2str(rngRxZeta(1)*180/pi,'%0.2f')];
hold all;
for idxZ = 2:length(rngRxZeta)
    plot(180*rngRxAlpha/pi,C_ZA(idxZ,:));
    lgd{idxZ} = ['\zeta = ' num2str(rngRxZeta(idxZ)*180/pi,'%0.2f')];
end
grid on;
xlabel('Azimuth Angle');
ylabel('Capacity (b/s/Hz)');
legend(lgd);
hold off;
axis auto;
if fSAVEALL
    fname = [ctDirRes 'fig_CapvsAzivsZ.fig'];
    saveas(gcf,fname,'mfig');
    fname = [ctDirRes 'CapvsAzivsZ.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% run('\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scratch\scratchSISO');

set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);