% scratch
close all;
clearvars;
clc;

% change responsivities of PDs !

%% SETUP
ctDirRes = 'C:\\Users\\pbutala\\Documents\\_Work\\Code V2\\Results\\2. l-MIMO\\20130521 400lx 4W\\';
ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\_Work\\Code V2\\Scratch\\scratchMIMOlambda.m';
ctFileCodeDest = [ctDirRes 'l-MIMOcode.m'];
ctFileVars = [ctDirRes 'l-MIMOdata.mat'];
fSAVEALL = false;

%% constants
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

Rpsd = getSOG(627,10,1,lambdas);
Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
Gpsd = getSOG(530,15,1,lambdas);
Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
Bpsd = getSOG(470,10,1,lambdas);
Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
Apsd = getSOG(590,10,1,lambdas);
Apsd = Apsd/(sum(Apsd)*LAMBDADELTA);

Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,4*Rpsd);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,4*Gpsd);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,4*Bpsd);
Ach = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,4*Apsd);

Wch = Rch + Gch + Bch + Ach;

Ambpsd = 5.8e-6*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

Rf = getSOG(627,10,1,lambdas);
Rf = 0.7*Rf/max(Rf);
Gf = getSOG(530,10,1,lambdas);
Gf = 0.7*Gf/max(Gf);
Bf = getSOG(470,10,1,lambdas);
Bf = 0.7*Bf/max(Bf);
Af = getSOG(590,10,1,lambdas);
Af = 0.7*Af/max(Af);

%% initialize room
room = cRoom;
[txX txY txZ] = meshgrid(room.L/2,room.W/2,3);
room.luminaire = cLuminaire([Rch Gch Bch Ach],txX,txY,txZ);

[rxX rxY rxZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
% [rxX rxY rxZ] = getGrid(1,1,1,1,1,1);
rxX = rxX+room.L/2;
rxY = rxY+room.W/2;
rxZ = rxZ + 1;
rxLoc = cLocation(rxX,rxY,rxZ);

Rrx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Rrx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
Grx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Grx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
Brx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Brx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
Arx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Arx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Af);

% for each receiver group
RH = room.getFreeSpaceGain(Rrx.location,Rrx.orientation,Rrx.rxFOV);
GH = room.getFreeSpaceGain(Grx.location,Grx.orientation,Grx.rxFOV);
BH = room.getFreeSpaceGain(Brx.location,Brx.orientation,Brx.rxFOV);
AH = room.getFreeSpaceGain(Arx.location,Arx.orientation,Arx.rxFOV);
% H(rx,lm,tx,ch);
RIsig = 0;GIsig = 0;BIsig = 0;AIsig = 0;
RIamb = 0;GIamb = 0;BIamb = 0;AIamb = 0;
nL = numel(room.luminaire);
for iL = 1:nL
    % TODO: calculate inter channel interference.
    color = room.luminaire(iL).color;
    channels = room.luminaire(iL).channels;
    % IR = rxR.getSignal(rxPsd);
%     [sig amb] = Rrx.getSignal(channels,sum(RH(:,iL,:,:),4)/4,Ambch);
    [sig amb] = Rrx.getSignal(channels,RH(:,iL,:,:),Ambch);
    RIsig = RIsig + sig;
    RIamb = RIamb + amb;
    % IG = rxG.getSignal(rxPsd);
%     [sig amb] = Grx.getSignal(color,sum(GH(:,iL,:,:),4)/4,Ambch);
    [sig amb] = Grx.getSignal(channels,GH(:,iL,:,:),Ambch);
    GIsig = GIsig + sig;
    GIamb = GIamb + amb;
    % IB = rxB.getSignal(rxPsd);
%     [sig amb] = Brx.getSignal(color,sum(BH(:,iL,:,:),4)/4,Ambch);
    [sig amb] = Brx.getSignal(channels,BH(:,iL,:,:),Ambch);
    BIsig = BIsig + sig;
    BIamb = BIamb + amb;
    % IA = rxA.getSignal(rxPsd);
%     [sig amb] = Arx.getSignal(color,sum(AH(:,iL,:,:),4)/4,Ambch);
    [sig amb] = Arx.getSignal(channels,AH(:,iL,:,:),Ambch);
    AIsig = AIsig + sig;
    AIamb = AIamb + amb;
end

%% RED channel
RIn = Rrx.tia.getNoise(40e6);
RIsh = 2*1.6e-19*(RIsig+RIamb)*40e6;
RSNR = (RIsig(:,:,:,1).^2)./(RIn+sum(RIsh,4));
RSNRdb = 10*log10(RSNR);

%% Green Channel
GIn = Grx.tia.getNoise(40e6);
GIsh = 2*1.6e-19*(GIsig+GIamb)*40e6;
% GSNR = (GIsig.^2)./(GIn+GIsh);
GSNR = (GIsig(:,:,:,2).^2)./(GIn+sum(GIsh,4));
GSNRdb = 10*log10(GSNR);

%% Blue Channel
BIn = Brx.tia.getNoise(40e6);
BIsh = 2*1.6e-19*(BIsig+BIamb)*40e6;
% BSNR = (BIsig.^2)./(BIn+BIsh);
BSNR = (BIsig(:,:,:,3).^2)./(BIn+sum(BIsh,4));
BSNRdb = 10*log10(BSNR);

%% Amber Channel
AIn = Arx.tia.getNoise(40e6);
AIsh = 2*1.6e-19*(AIsig+AIamb)*40e6;
% ASNR = (AIsig.^2)./(AIn+AIsh);
ASNR = (AIsig(:,:,:,4).^2)./(AIn+sum(AIsh,4));
ASNRdb = 10*log10(ASNR);

%% FIGURES
if fSAVEALL
    delete([ctDirRes '*.png']);
end
ifig = 1;
% draw setup
figure;
room.drawSetup(cLocation(room.L/2,room.W/2,rxZ(1)),cOrientation(0,0,0),Rrx.rxFOV);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Setup.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

% plot LED PSDs
figure;
plot(Rch.npsd.X,Rch.npsd.Y/Rch.npsd.Ymax,'r');
hold on;
plot(Gch.npsd.X,Gch.npsd.Y/Gch.npsd.Ymax,'g');
plot(Bch.npsd.X,Bch.npsd.Y/Bch.npsd.Ymax,'b');
plot(Ach.npsd.X,Ach.npsd.Y/Ach.npsd.Ymax,'y');
grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized PSD');
% title('Normalized PSDs of LEDs');
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' LED PSD.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

figure;
plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax,'k');
xlabel('Wavelength (nm)');
ylabel('Normalized PSD');
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Total PSD.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

% CIE plot
figure;
obs = cCIE;
obs.getCoordinates([Rch.npsd Gch.npsd Bch.npsd Ach.npsd Wch.npsd]);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' CIEplot.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

% plot filter responses
figure;
plot(Rrx.sensor.filter.X,Rrx.sensor.filter.Y,'r');
hold on;
plot(Grx.sensor.filter.X,Grx.sensor.filter.Y,'g');
plot(Brx.sensor.filter.X,Brx.sensor.filter.Y,'b');
plot(Arx.sensor.filter.X,Arx.sensor.filter.Y,'y');
title('Filter responses');
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' FilterResponse.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

% plot illuminance,irradiance
figure;
room.drawIlluminance;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Illuminance.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

figure;
room.drawIrradiance;
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' Irradiance.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

cLim = [min([RSNRdb(:);GSNRdb(:);BSNRdb(:);ASNRdb(:)]) ...
    max([RSNRdb(:);GSNRdb(:);BSNRdb(:);ASNRdb(:)])];

% plot RGBA SNRs
figure;
subplot(2,2,1);
set(gca,'climmode','manual');
set(gca,'clim',cLim);
surf(rxX,rxY,reshape(RSNRdb,size(rxX)));
axis([0 room.L 0 room.W min(RSNRdb(:)) max(RSNRdb(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('SNR (Red)\nMax: %0.2f dB',max(RSNRdb(:)));
title(tStr);

subplot(2,2,2);
set(gca,'climmode','manual');
set(gca,'CLim',cLim);
surf(rxX,rxY,reshape(GSNRdb,size(rxX)));
axis([0 room.L 0 room.W min(GSNRdb(:)) max(GSNRdb(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('SNR (Green)\nMax: %0.2f dB',max(GSNRdb(:)));
title(tStr);

subplot(2,2,3);
set(gca,'climmode','manual');
set(gca,'CLim',cLim);
surf(rxX,rxY,reshape(BSNRdb,size(rxX)));
axis([0 room.L 0 room.W min(BSNRdb(:)) max(BSNRdb(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('SNR (Blue)\nMax: %0.2f dB',max(BSNRdb(:)));
title(tStr);

subplot(2,2,4);
set(gca,'climmode','manual');
set(gca,'CLim',cLim);
surf(rxX,rxY,reshape(ASNRdb,size(rxX)));
axis([0 room.L 0 room.W min(ASNRdb(:)) max(ASNRdb(:))]);
grid on;
colorbar;
view(3);
tStr = sprintf('SNR (Amber)\nMax: %0.2f dB',max(ASNRdb(:)));
title(tStr);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' SNR.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end

if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest);
    save(ctFileVars);
end
%% code for colored reflections
% for iR = 1:rxR.rxCount;
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
%     ambPSD(iR,1) = chAmb;
% end
% [IRsig IRamb] = rxR.getSignal(rxPSD,ambPSD);


























