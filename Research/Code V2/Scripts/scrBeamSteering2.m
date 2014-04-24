% scrBeamSteering
% close all;
scrBeamSteering;
% figure;
clearvars -except figSNR ITER;
clc;

% DEFAULT COSMETIC SETTINGS
dlinelw = get(0,'DefaultLineLineWidth');
set(0,'DefaultLineLineWidth',2);
daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',10);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');
dfigppm = get(0,'DefaultFigurePaperPositionMode');
set(0,'DefaultFigurePaperPositionMode','Manual');
dfigpu = get(0,'DefaultFigurePaperUnits');
set(0,'DefaultFigurePaperUnits','Inches');
dfigpp = get(0,'DefaultFigurePaperPosition');
set(0,'DefaultFigurePaperPosition',[0 0 11 8.5]);

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fFILTBLUE = true;
fNOISESHOTONLY = false;
fARCHIVE = true;

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = '_sc_m20_m20_steer';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%%  Set variables
varIllSet = 400; % varILLSet
varLkD = 2;
varTxa = 1e-2;
varRxa = 5e-3;
varRxZeta = 0;
varRxAlpha = 0;
varRxTau = 0;
varLm = 20; % m=20 => ~15 deg half angle
% varLm = -log(2)/log(cos(pi/12));
rdFlMax = realmin('double');
rdFlMin = realmax('double');
IFIG = 1;
figILL = [];
figure(figSNR);
subplot(3,1,[1 2]);
[~,~,LGDH,LGD] = legend;

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\BeamSteering2\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrBeamSteering2.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/BeamSteering2/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrBeamSteering2.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\MatlabResults\\BeamSteering2\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrBeamSteering2.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrBeamSteering2' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'scrBeamSteering2' CHARIDXARCHIVE '.mat'];
ctrxB = 40e6;
ctILLTh = 300;

%% constants
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
% etc vars
clLC = {'k','b','r','g','m'};
lenLC = length(clLC);
clLS = {'-.','--',':'};
lenLS = length(clLS);
clMK = {'h','o','s','d','v','^','<','>','p'};
lenMK = length(clMK);

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
Wch1 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
Wch2 = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);

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
txX = (-0.75:0.5:0.75) + room.L/2;
txY = zeros(size(txX)) + room.W/2;
txZ = zeros(size(txX)) + 3;
txLoc = cLocation(txX,txY,txZ);
room.luminaire = cLuminaire(Wch,txLoc);
room.luminaire.order = 1;
room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
LMFLUX = room.luminaire.lmFlux;

rxX = 0:0.1:room.L;
rxY = zeros(size(rxX))+ room.W/2;
rxZ = zeros(size(rxX)) + 3 - varLkD;

[ilX,ilY,ilZ] = getGrid(room.L,room.W,1,0.1,0.1,2,'Fill');
ilX = ilX+room.L/2;
ilY = ilY+room.W/2;
ilZ = ilZ + 3 - varLkD;

SNRm = zeros(numel(rxX),numel(txX),numel(txX));
% for txj=1:room.lmCount
for txj=1:numel(txX)
    for rxi=1:numel(rxX)
        clear room lm1 lm2;
        room = cRoom(4,4,4);
        txLoc1 = cLocation(txX(1:end ~=txj),txY(1:end ~=txj),txZ(1:end ~=txj));
        txLoc2 = cLocation(txX(txj),txY(txj),txZ(txj));
        lm1 = cLuminaire(Wch1,txLoc1);
        lm1.order = 20;
        lm2 = cLuminaire(Wch2,txLoc2);
        lm2.order = varLm;
        room.luminaire = [lm1 lm2];
%         room.luminaire = [cLuminaire(Wch,txLoc1) cLuminaire(Wch,txLoc2)];
%         room.luminaire(1).scaleOutputFlux(LMFLUX/room.luminaire(1).lmFlux);
%         room.luminaire(2).scaleOutputFlux(LMFLUX/room.luminaire(2).lmFlux);
%         room.luminaire(2).order = varLm;
        
        rcvr = cSinglePixelReceiverWhiteReflection(rxX(rxi),rxY(rxi),rxZ(rxi));
        rxOri = cOrientation(varRxZeta,varRxAlpha,varRxTau);
        
        if fFILTBLUE
            rcvr.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
            rcvr.sensor.responsivity = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
        end
        
        rcvr.orientation = rxOri;
        rcvr.sensor.dimension.L = varRxa;
        rcvr.sensor.dimension.W = varRxa;
        rcvr.sensor.dimension.H = 1e-3;
        rcvr.optics.fov = pi/2;
        % beamsteer
        Vect = [txX(txj)-rxX(rxi);txY(txj)-rxY(rxi);txZ(txj)-rxZ(rxi)];
        VLen = sqrt(sum(Vect.^2));
        Z = pi-acos(Vect(3)/VLen);
        A = atan2(-Vect(1),-Vect(2));
        txOri2 = cOrientation(Z,A,0);
        room.luminaire(2).orientation = txOri2;
        
        % set illumination
%         room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
        room.setIlluminance(cLocation(ilX,ilY,ilZ),varIllSet);
        if room.luminaire(1).rdFlux > rdFlMax
            rdFlMax = room.luminaire(1).rdFlux;
        end
        if room.luminaire(1).rdFlux < rdFlMin
            rdFlMin = room.luminaire(1).rdFlux;
        end
        % for each receiver group
        H = room.getFreeSpaceGain(rcvr.location,rcvr.orientation,rcvr.rxFOV);
        nL = numel(room.luminaire);
        Isig = zeros(1,sum(room.lmCount));
        
        color = room.luminaire(1).color;
        [sig1,Iamb] = rcvr.getSignal(color,sum(H(:,1,:),4),Ambch);
        Isig(1:end ~=txj) = sig1;
        
        color = room.luminaire(2).color;
        [sig2,~] = rcvr.getSignal(color,sum(H(:,2,:),4),Ambch);
        Isig(txj) = sig2(1);
        
        if fNOISESHOTONLY
            In = 0;
        else
            In = rcvr.tia.getNoise(ctrxB);
        end
        
        Ish = 2*1.6e-19*(Isig+Iamb)*ctrxB;
        for j=1:numel(txX)
            SNRm(rxi,j,txj) = (Isig(j)^2)/(sum(Ish(:))+In);
        end
%         if (txj==1) && (rxi==ceil(numel(rxX)/2))
        if (txj==1) && (rxi==13)
            figure(figSNR);
            subplot(3,3,9);
            room.drawIlluminance(cLocation(ilX,ilY,ilZ),ctILLTh,gca);    
            figure(figSNR);
            subplot(3,3,9);
            caxis([0 2500]);
            axis([0 room.L 0 room.W 0 2500]);
            title(sprintf('Illumination (steer)\ntx:[%0.2f %0.1f %0.1f] rx:[%0.1f %0.1f %0.1f]',...
                txX(txj),txY(txj),txZ(txj),rxX(rxi),rxY(rxi),rxZ(rxi)));
            
            room.drawIlluminance(cLocation(ilX,ilY,ilZ));
            figILL(end+1) = gcf;
            caxis([0 2500]);
            axis([0 room.L 0 room.W 0 2500]);
            title(sprintf('Illumination (steer)\ntx:[%0.2f %0.1f %0.1f] rx:[%0.1f %0.1f %0.1f]',...
                txX(txj),txY(txj),txZ(txj),rxX(rxi),rxY(rxi),rxZ(rxi)));
            
            % draw setup
            figSetup = figure;
            room.drawSetup(cLocation(rxX(rxi),rxY(rxi),rxZ(rxi)),cOrientation(0,0,0),pi/4);
            if fSAVEALL
                fname = [ctDirRes num2str(IFIG) ' Setup' CHARIDXARCHIVE];
                saveas(gcf,[fname '.fig'],'fig');
                saveas(gcf,[fname '.png'],'png');
                IFIG = IFIG + 1;
            end
            rotate3d on;
        end
    end
end
SNRmdb = 10*log10(SNRm);
%% FIGURES
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end

% figSNR = figure;
figure(figSNR);
subplot(3,1,[1 2]);
% hold all;
iLS = rem(ITER,lenLS)+1;
% for j=1:numel(txX)
%     iLC = rem(j,lenLC)+1;
%     iMK = rem(j,lenMK)+1;
%     plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%     plot(rxX,max(SNRmdb(:,:,j),[],2),plStyle);
% end
% iMK = rem(j+1,lenMK)+1;
iLC = rem(ITER,lenLC)+1;
iMK = rem(ITER,lenMK)+1;
% plStyle = ['k' clLS{iLS} clMK{iMK}];
plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
LGDH(end+1) = plot(rxX,max(reshape(SNRmdb,numel(rxX),numel(txX)*numel(txX)),[],2),plStyle);
LGD{end+1} = sprintf('steer %-2.2f:%-2.2f W',rdFlMin,rdFlMax);
% axis([min(rxX) max(rxX) 25 2*ceil(max(SNRmdb(:))/2)]);
axis([min(rxX) max(rxX) 25 50]);

legend(LGDH,LGD);
grid on;
if fSAVEALL
    fname = [ctDirRes num2str(IFIG) ' SNRbmstr' CHARIDXARCHIVE];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    IFIG = IFIG + 1;
    
    for f=1:numel(figILL)
        fname = [ctDirRes num2str(IFIG) ' illumination' CHARIDXARCHIVE];
        saveas(figILL(f),[fname '.fig'],'fig');
        saveas(figILL(f),[fname '.png'],'png');
        IFIG = IFIG + 1;
    end
end

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
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);





















