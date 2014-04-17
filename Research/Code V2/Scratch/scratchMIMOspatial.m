% scratch
close all;
clearvars;
clc;

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = true;
fFILTBLUE = false;
fNOISESHOTONLY = true;

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\3. s-MIMO\0. Scratch\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scratch\scratchMIMOspatial.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/3. s-MIMO/0. Scratch/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scratch/scratchMIMOspatial.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\_Work\\Code V2\\Results\\3. s-MIMO\\0. Scratch\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\_Work\\Code V2\\Scratch\\scratchMIMOspatial.m';
    otherwise
        error('Station not defined');
end

ctFileCodeDest = [ctDirRes 's-MIMO.m'];
ctFileVars = [ctDirRes 's-MIMOdata.mat'];

%% constants
ctILLTh = 300;
ctILLSet = 400;
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
% room.setIlluminance(.) will set the output power of the transmitters

Ambpsd = 5.8e-6*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

if fFILTBLUE
    Bf = getSOG(470,10,1,lambdas);
    Bf = 0.7*Bf/max(Bf);
end

f = 5e-3;
% f = 1.38e-2;
% fN = f/sqrt(4e-6/pi);
fN = 1;
% M = f/(2-f);
% D = 5e-3/M;

%% initialize room
room = cRoom(4,4,4);
% [txX txY txZ] = getGrid(2,2,1,0.5,0.5,2);
[txX txY txZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
txX = txX + room.L/2;
txY = txY + room.W/2;
txZ = txZ + 3;
txLoc = cLocation(txX,txY,txZ);
txOri = cOrientation(pi,0,0);
room.luminaire = cLuminaire(Wch,txLoc);
room.luminaire.orientation = txOri;
% room.luminaire.dimension = cSize(36e-2,36e-2,1e-3);
room.setIlluminance(cLocation(room.L/2,room.W/2,1),ctILLSet);

% [rxX rxY rxZ] = getGrid(room.L,room.W,1,1,1,2,'Fill');
[rxX rxY rxZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
rxX = rxX + room.L/2;
rxY = rxY + room.W/2;
rxZ = rxZ + 1;
% [rxX rxY rxZ] = getGrid(1,1,1);
rxLoc = cLocation(rxX,rxY,rxZ);

NPxX = 10;
NPxY = 10;
NPxD = 1e-3;
% NPxD = 5.91e-2/NPxX;
[pxX pxY pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);
Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,f,fN);
Wrx.optics.T = 0.9; % set lens transmission to 90%
Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;

if fFILTBLUE
    Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
    Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
end

% for each receiver group
H = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);
nL = numel(room.luminaire);
txSz = room.lmDimension;

% currently assuming just 1 lum group at multiple loc
color = room.luminaire(1).color;
txOri = room.luminaire(1).orientation;
[Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = ...
    Wrx.getSignal(color,squeeze(sum(H(:,1,:,:),4)/4),Ambch,txLoc,txSz,txOri); 

if fNOISESHOTONLY
    In = 0;
else
    In = Wrx.tia.getNoise(40e6);
end
Ish = 2*1.6e-19*(Isig+Iamb)*40e6;
SNR = (Isig.^2)./(In(1)+Ish);
SNRdb = 10*log10(SNR);
cLim = [floor(min(SNRdb(:))) ceil(max(SNRdb(:)))];

%% SAVE workspace
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest);
    save(ctFileVars);
end

%% FIGURES
if fSAVEALL
    delete([ctDirRes '*.png']);
end
ifig = 1;

% draw setup
locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
figure;
room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Setup.png'];
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
    fname = [ctDirRes num2str(ifig) ' LED PSD.png'];
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
    fname = [ctDirRes num2str(ifig) ' CIEplot.png'];
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
    fname = [ctDirRes num2str(ifig) ' FilterResponse.png'];
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
    fname = [ctDirRes num2str(ifig) ' SensorResponsivity.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% draw spots
% loc2 = cLocation(room.L/2+0.25,room.W/2,rxZ(1));
figure;
h1 = gca;
% h2 = subplot(1,2,2);
% Wrx.getImage(locCntr+loc2,txLoc,txSz,txOri,[h1 h2]);
Wrx.getImage(locCntr,txLoc,txSz,txOri,h1);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Spots.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% plot illuminance,irradiance
% [plX plY plZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
% plLoc = cLocation(plX+room.L/2,plY+room.W/2,plZ+0.85);
figure;
H(1) = gca;
rotate3d on;
figure;
H(2) = gca;
rotate3d on;
room.drawIlluminance(Wrx.location,ctILLTh,H);
if fSAVEALL
    figure(get(H(1),'Parent'));
    fname = [ctDirRes num2str(ifig) ' Illuminance 3D.png'];
    saveas(gcf,fname,'png');
    view(2);
    fname = [ctDirRes num2str(ifig) ' Illuminance 2D.png'];
    saveas(gcf,fname,'png');
    view(3);
    if fCLOSEFIGS
        close;
    end
    ifig = ifig + 1;
    figure(get(H(2),'Parent'));
    fname = [ctDirRes num2str(ifig) ' Illuminance Coverage.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
    if fCLOSEFIGS
        close;
    end
end
figure;
room.drawIrradiance;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Irradiance 3D.png'];
    saveas(gcf,fname,'png');
    view(2);
    fname = [ctDirRes num2str(ifig) ' Irradiance 2D.png'];
    saveas(gcf,fname,'png');
    view(3);
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
else
    rotate3d on;
end

% plot SNRS
% for each transmitter calculate SNR mesh and plot it
[nRx nTx] = size(SNRdb);
for iT = 1:nTx
    snrdb = SNRdb(:,iT);
    figure;
    set(gca,'climmode','manual');
    set(gca,'CLim',cLim);
    surf(rxX,rxY,reshape(snrdb,size(rxX)),'FaceColor','interp');
    axis([0 room.L 0 room.W min(snrdb(:)) max(snrdb(:))]);
    grid on;
    colorbar;
    view(3);
    tStr = sprintf('SNR for Tx @ (%0.2f %0.2f %0.2f)\nMax: %0.2f dB',txLoc.X(iT),txLoc.Y(iT),txLoc.Z(iT),max(snrdb(:)));
    title(tStr);
    xlabel('X');
    ylabel('Y');
    if fSAVEALL
        fname = [ctDirRes num2str(ifig) ' SNR Tx ' num2str(txLoc.X(iT),'%0.2f') ' ' num2str(txLoc.Y(iT),'%0.2f') ' ' num2str(txLoc.Z(iT),'%0.2f') ' 3D.png'];
        saveas(gcf,fname,'png');
        view(2);
        fname = [ctDirRes num2str(ifig) ' SNR Tx ' num2str(txLoc.X(iT),'%0.2f') ' ' num2str(txLoc.Y(iT),'%0.2f') ' ' num2str(txLoc.Z(iT),'%0.2f') ' 2D.png'];
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

% other

%% colored code
% for iR = 1:Wrx.rxCount;
%     iTx = 1;
%     for iL = 1:nL % lumianire group
%         ch = room.luminaire(iL).channels;
%         % for each transmitter in the luminaire group
%         for iT = 1:room.luminaire(iL).lmCount;
%             h = squeeze(H(iR,iL,iT,:));
%             rxPSD(iR,iTx) = ch*h;
%             iTx = iTx + 1;
%         end
%     end
%     ambPSD(iR) = Ambch;
% end
    
% [Ipx Iambpx Isig Iamb frSpOnPx frPxOnSp] = Wrx.getSignal(rxPSD,ambPSD,txLoc,txSz); % currently assuming just 1 lum group at multiple loc

% In = Wrx.tia.getNoise(40e6);
% Ish = 2*1.6e-19*(Isig+Iamb)*40e6;
% 
% SNR = (Isig.^2)./(In(1)+Ish);
% SNRdb = 10*log10(SNR);

% [txM txN] = size(txLoc.X);
% IsigMt = zeros([txM txN rxW.rxCount]);
% for id = 1:rxW.rxCount
%     IsigMt(1:txM,1:txN,id) = reshape(Isig(id,:),txM,txN);
% end
% IpxMt = zeros(NPxX,NPxY,rxW.rxCount);
% for id = 1:rxW.rxCount
%     IpxMt(1:NPxX,1:NPxY,id) = reshape(Ipx(id,:),NPxX,NPxY);
% end
% IambMt = zeros([txM txN rxW.rxCount]);
% for id = 1:rxW.rxCount
%     IambMt(1:txM,1:txN,id) = reshape(Iamb(id,:),txM,txN);
% end
% IambpxMt = zeros(NPxX,NPxY,rxW.rxCount);
% for id = 1:rxW.rxCount
%     IambpxMt(1:NPxX,1:NPxY,id) = reshape(Iambpx(id,:),NPxX,NPxY);
% end
% display(IsigMt);
% display(IpxMt);
% display(IambMt);
% display(IambpxMt);

% SNR = IsigMt^2/(In+Ish);
% SNRdb = 10*log10(SNR);

% for iT = 1:size(SNRdb,2) % for each transmitter
%     rxSNRdb = reshape(SNRdb(:,iT),size(rxX));
%     figure;
%     surf(rxX,rxY,rxSNRdb);
%     axis([0 room.L 0 room.W min(rxSNRdb(:)) max(rxSNRdb(:))]);
%     grid on;
%     colorbar;
%     view(3);
%     tStr = sprintf('SNR tx@(%0.2f,%0.2f,%0.2f)\nMax: %0.2f dB',txX(iT),txY(iT),txZ(iT),max(rxSNRdb(:)));
%     title(tStr);
% end
























