% script SISO
close all;
clearvars;
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');

%% FLAGS
fSTATION = 3;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fFILTBLUE = false;
fNOISESHOTONLY = false;

%% VARIABLE PARAMS

%%  Set variables
varIllSet = 400; % varILLSet
varLkD = 2;
varTxa = 1e-2;
varRxa = 5e-3;
varRxZeta = 0;
varRxAlpha = 0;
varRxTau = 0;

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\1. SISO\1. Basic\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrSISO.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/1. SISO/1. Basic/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrSISO.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\1. SISO\\1. Basic\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrSISO.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrSISO.m'];
ctFileVars = [ctDirRes 'SISOdata.mat'];
ctrxB = 40e6;
ctILLTh = 200;

%% constants
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);

Ambpsd = 5.8e-6*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

if fFILTBLUE
    % Bf = getSOG(470,10,1,lambdas);
    % Bf = 0.7*Bf/max(Bf);
    Bf = zeros(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
end

%% initialize room
room = cRoom(4,4,4);
txLoc = cLocation(room.L/2,room.W/2,3);
room.luminaire = cLuminaire(Wch,txLoc);
room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);

[rxX,rxY,rxZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');

rxX = rxX+room.L/2;
rxY = rxY+room.W/2;
rxZ = rxZ + 3 - varLkD;
rcvr = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
rxOri = cOrientation(varRxZeta,varRxAlpha,varRxTau);

if fFILTBLUE
    rcvr.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
    rcvr.sensor.responsivity = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
end

rcvr.orientation = rxOri;
rcvr.sensor.dimension.L = varRxa;
rcvr.sensor.dimension.W = varRxa;
rcvr.sensor.dimension.H = 1e-3;
rcvr.optics.fov = pi/3;

% for each receiver group
H = room.getFreeSpaceGain(rcvr.location,rcvr.orientation,rcvr.rxFOV);
nL = numel(room.luminaire);
Isig = 0;
Iamb = 0;
for iL = 1:nL
    color = room.luminaire(iL).color;
    [sig,amb] = rcvr.getSignal(color,sum(H(:,iL,:),4),Ambch);
    Isig = Isig + sig;
    Iamb = Iamb + amb;
end

if fNOISESHOTONLY
    In = 0;
else
    In = rcvr.tia.getNoise(ctrxB);
end

Ish = 2*1.6e-19*(Isig+Iamb)*ctrxB;
SNR = (Isig.^2)./(In+Ish);
SNRdb = 10*log10(SNR);
C_siso = log2(1 + SNR);

%% FIGURES
if fSAVEALL
    delete([ctDirRes '*.png']);
    delete([ctDirRes '*.mat']);
    delete([ctDirRes '*.fig']);
end
ifig = 1;

% draw setup
figSetup = figure;
room.drawSetup(cLocation(room.L/2,room.W/2,rxZ(1)),cOrientation(0,0,0),rcvr.rxFOV);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Setup'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end
rotate3d on;

% plot LED PSDs
figPSD = figure;
plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax);
grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized PSD');
% title('Normalized PSDs of LEDs');
axis tight;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' LED PSD'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end

% CIE plot
figCIE = figure;
obs = cCIE;
obs.getCoordinates(Wch.npsd);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' CIEplot'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end

% plot filter responses
figFilt = figure;
plot(rcvr.sensor.filter(1).X,rcvr.sensor.filter(1).Y);
% title('Filter Transmission');
grid on;
xlabel('Wavelength (nm)');
ylabel('Filter Transmission');
axis tight;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' FilterResponse'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end

% plot receiver responses
figResp = figure;
plot(rcvr.sensor.responsivity(1).X,rcvr.sensor.responsivity(1).Y);
% title('Sensor Responsivity');
grid on;
xlabel('Wavelength (nm)');
ylabel('Responsivity (A.W^{-1})');
axis tight;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' SensorResponsivity'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end

% plot illuminance,irradiance
figIll_1 = figure;
H(1) = gca;
rotate3d on;
figIll_2 = figure;
H(2) = gca;
rotate3d on;
room.drawIlluminance(rcvr.location,ctILLTh,H);
if fSAVEALL
    figure(get(H(1),'Parent'));
    fname = [ctDirRes num2str(ifig) ' Illuminance'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname ' 3D.png'],'png');
    view(2);
    saveas(gcf,[fname ' 2D.png'],'png');
    view(3);
    ifig = ifig + 1;

    figure(get(H(2),'Parent'));
    fname = [ctDirRes num2str(ifig) ' Illuminance Coverage'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
end
figIrr = figure;
room.drawIrradiance;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Irradiance'];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname ' 3D.png'],'png');
    view(2);
    saveas(gcf,[fname ' 2D.png'],'png');
    view(3);
    ifig = ifig + 1;
end
rotate3d on;

% plot SNR
figSNR = figure;
surf(rxX,rxY,reshape(SNRdb,size(rxX)),'FaceColor','interp');
axis([0 room.L 0 room.W min(SNRdb(:)) max(SNRdb(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('Max: %0.2f dB',max(SNRdb(:)));
title(tStr);
xlabel('X');
ylabel('Y');
zStr = sprintf('SNR (dB)');
zlabel(zStr);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' SNR Tx ' num2str(txLoc.X,'%0.2f') ' ' num2str(txLoc.Y,'%0.2f') ' ' num2str(txLoc.Z,'%0.2f')];
    saveas(gcf,[fname ' .fig'],'fig');
    saveas(gcf,[fname ' 3D.png'],'png');
    view(2);
    saveas(gcf,[fname ' 2D.png'],'png');
    view(3);
    ifig = ifig + 1;
end
rotate3d on;

% plot Capacity
figCAP = figure;
surf(rxX,rxY,reshape(C_siso,size(rxX)),'FaceColor','interp');
axis([0 room.L 0 room.W min(C_siso(:)) max(C_siso(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('Max: %0.2f b.s^{-1}.Hz^{-1}',max(C_siso(:)));
title(tStr);
xlabel('X');
ylabel('Y');
zStr = sprintf('Capacity (b.s^{-1}.Hz^{-1})');
zlabel(zStr);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' CAP Tx ' num2str(txLoc.X,'%0.2f') ' ' num2str(txLoc.Y,'%0.2f') ' ' num2str(txLoc.Z,'%0.2f')];
    saveas(gcf,[fname ' .fig'],'fig');
    saveas(gcf,[fname ' 3D.png'],'png');
    view(2);
    saveas(gcf,[fname ' 2D.png'],'png');
    view(3);
    ifig = ifig + 1;
end
rotate3d on;
% other
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest);
    save(ctFileVars);
end

% restore defaults
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);






















