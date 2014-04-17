% scratch
close all;
clearvars;
clc;

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

Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);
Ach = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Apsd);
chTotal = Rch + Gch + Bch + Ach;

Rf = getSOG(627,10,1,lambdas);
Rf = 0.7*Rf/max(Rf);
Gf = getSOG(530,10,1,lambdas);
Gf = 0.7*Gf/max(Gf);
Bf = getSOG(470,10,1,lambdas);
Bf = 0.7*Bf/max(Bf);
Af = getSOG(590,10,1,lambdas);
Af = 0.7*Af/max(Af);

Ambpsd = 5.8e-6*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

f = 3e-2;
fN = f/sqrt(4e-6/pi);
M = f/(2-f);
% D = 5e-3/M;

%% initialize room
room = cRoom(7,4,4);
[txX txY txZ] = getGrid(3,3,1,0.5,0.5,2);
% [txX txY txZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
txLoc = cLocation(txX+room.L/2,txY+room.W/2,txZ+3);
txOri = cOrientation(pi,0,0);
room.luminaire = cLuminaire([Rch Gch Bch Ach],txLoc);
room.luminaire.orientation = txOri;

% [rxX rxY rxZ] = getGrid(room.L,room.W,1,2,2,2,'Fill');
[rxX rxY rxZ] = getGrid(1,1,1);
rxLoc = cLocation(rxX+room.L/2,rxY+room.W/2,rxZ+1);

NPxX = 20;
NPxY = 20;
NPxD = 1e-3;
[pxX pxY pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);

rcvr = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,f,fN);

% SET SENSOR FILTERS

rcvr.sensor.filter(1:2:NPxX,1:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
rcvr.sensor.filter(1:2:NPxX,2:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
rcvr.sensor.filter(2:2:NPxX,1:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
rcvr.sensor.filter(2:2:NPxX,2:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Af);
rcvr.sensor.type(1:2:end,1:2:end) = 1;
rcvr.sensor.type(1:2:end,2:2:end) = 2;
rcvr.sensor.type(2:2:end,1:2:end) = 3;
rcvr.sensor.type(2:2:end,2:2:end) = 4;
rcvr.sensor.color = [1 0 0;0 1 0;0 0 1;1 1 0];
% for each receiver group
H = room.getFreeSpaceGain(rcvr.location,rcvr.orientation,rcvr.rxFOV);
nL = numel(room.luminaire);

txSz = room.lmDimension;
color = room.luminaire(1).color;
txOri = room.luminaire(1).orientation;
[Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = rcvr.getSignal(color,sum(H(:,1,:,:),4)/4,Ambch,txLoc,txSz,txOri); % currently assuming just 1 lum group at multiple loc

In = rcvr.tia.getNoise(40e6);
Ish = 2*1.6e-19*(Isig+Iamb)*40e6;
SNR = (Isig.^2)./(In(1)+Ish);
SNRdb = 10*log10(SNR);

%% FIGURES
locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
figure;
room.drawSetup(locCntr,cOrientation(0,0,0),rcvr.rxFOV);

h = figure;
h1 = gca;
locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
rcvr.getImage(locCntr,txLoc,txSz,txOri,h1);

[plX plY plZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
plLoc = cLocation(plX+room.L/2,plY+room.W/2,plZ+0.85);
figure;
room.drawIlluminance(plLoc);
figure;
room.drawIrradiance(plLoc);


%% Colored code
% rcvr.sensor.filter(1:2:NPxX,1:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
% rcvr.sensor.filter(1:2:NPxX,2:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
% rcvr.sensor.filter(2:2:NPxX,1:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
% rcvr.sensor.filter(2:2:NPxX,2:2:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Af);
% rcvr.sensor.type(1:2:end,1:2:end) = 1;
% rcvr.sensor.type(1:2:end,2:2:end) = 2;
% rcvr.sensor.type(2:2:end,1:2:end) = 3;
% rcvr.sensor.type(2:2:end,2:2:end) = 4;
% rcvr = cImagingReceiverColoredReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,f,fN);
% for iR = 1:rcvr.rxCount;
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
    
% [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = rcvr.getSignal(rxPSD,ambPSD,txLoc,txSz); % currently assuming just 1 lum group at multiple loc
% ASSUMING JUST 1 LUMINAIRE TYPE

% Iamb = rxW.getSignal(Ambch,ambLoc,ambSz);
% In = rcvr.tia.getNoise(40e6);
% Ish = 2*1.6e-19*(Isig+Iamb)*40e6;
% Ish = 2*1.6e-19*Is*40e6;
% nR = rcvr.rxCount;
% nP = rcvr.sensor.location.count;
% nT = numel(txLoc.X);
% nC = numel(room.luminaire.channels);
% [txM txN] = size(txLoc.X);
% IsigMt = zeros([txM txN 4 rcvr.rxCount]);
% IsigWMt = zeros([txM txN 4 rcvr.rxCount]);
% for id = 1:rcvr.rxCount
%     IsigMt(1:txM,1:txN,:,id) = reshape(Isig(id,:),txM,txN,4);
%     IsigWMt(1:txM,1:txN,:,id) = reshape(IsigW(id,:),txM,txN,4);
% end
% 
% IpxMt = zeros(NPxX,NPxY,rcvr.rxCount);
% IpxWMt = zeros(NPxX,NPxY,rcvr.rxCount);
% for id = 1:rcvr.rxCount
%     IpxMt(1:NPxX,1:NPxY,id) = reshape(Ipx(id,:),NPxX,NPxY);
%     IpxWMt(1:NPxX,1:NPxY,id) = reshape(IpxW(id,:),NPxX,NPxY);
% end
% 
% IambMt = zeros([txM txN 4 rcvr.rxCount]);
% IambWMt = zeros([txM txN 4 rcvr.rxCount]);
% for id = 1:rcvr.rxCount
%     IambMt(1:txM,1:txN,:,id) = reshape(Iamb(id,:),txM,txN,4);
%     IambWMt(1:txM,1:txN,:,id) = reshape(IambW(id,:),txM,txN,4);
% end
% 
% IambpxMt = zeros(NPxX,NPxY,rcvr.rxCount);
% IambpxWMt = zeros(NPxX,NPxY,rcvr.rxCount);
% for id = 1:rcvr.rxCount
%     IambpxMt(1:NPxX,1:NPxY,id) = reshape(IambpxW(id,:),NPxX,NPxY);
%     IambpxWMt(1:NPxX,1:NPxY,id) = reshape(IambpxW(id,:),NPxX,NPxY);
% end
% display(IsigMt);
% display(IsigWMt);
% display(IpxMt);
% display(IpxWMt);
% display(IambMt);
% display(IambpxMt);

% 

























