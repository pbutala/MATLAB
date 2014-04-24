% scrBeamSteering
close all;
clearvars;
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
fARCHIVE = false;

CHAROVERWRITE = '~';
if(fARCHIVE)
    CHARIDXARCHIVE = '_sc_m1_m20_400lx';           % ARCHIVE INDEX
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
varLm = [1 20 20]; % m=20 => ~15 deg half angle
fIll400 = [true false true];
fPLOT = [true false true];
LGD = [];
LGDH = [];
ITER = 1;
%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\BeamSteering\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrBeamSteering.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/BeamSteering/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrBeamSteering.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\MatlabResults\\BeamSteering\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrBeamSteering.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrBeamSteering' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes 'scrBeamSteering' CHARIDXARCHIVE '.mat'];
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
room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
LMFLUX = room.luminaire.lmFlux;

rxX = 0:0.1:room.L;
rxY = zeros(size(rxX))+ room.W/2;
rxZ = zeros(size(rxX)) + 3 - varLkD;
rcvr = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
rxOri = cOrientation(varRxZeta,varRxAlpha,varRxTau);

[ilX,ilY,ilZ] = getGrid(room.L,room.W,1,0.1,0.1,2,'Fill');
ilX = ilX+room.L/2;
ilY = ilY+room.W/2;
ilZ = ilZ + 3 - varLkD;
if fFILTBLUE
    rcvr.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
    rcvr.sensor.responsivity = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
end

rcvr.orientation = rxOri;
rcvr.sensor.dimension.L = varRxa;
rcvr.sensor.dimension.W = varRxa;
rcvr.sensor.dimension.H = 1e-3;
rcvr.optics.fov = pi/2;

%% FIGURES
if fSAVEALL
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
end
ifig = 1;
iSF = 1;
figSNR = figure;
figILL = [];
subplot(3,1,[1 2]);
ylabel('SNR(db)');
xlabel(sprintf('Location X (Y=%d, Z=%d)',room.W/2,3 - varLkD));
hold on;
for idm = 1:numel(varLm)
    txm = varLm(idm);
    %% initialize room
    room.luminaire.order = txm;
    if fIll400(idm)
        room.setIlluminance(cLocation(ilX,ilY,ilZ),varIllSet);
    end
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
    
    SNR = zeros(rcvr.rxCount,room.lmCount);
    for j=1:room.lmCount
        for i=1:rcvr.rxCount
            SNR(i,j) = (Isig(i,1,j)^2)/(sum(Ish(i,1,:),3)+In);
        end
    end
    SNRdb = 10*log10(SNR);
    figure(figSNR);
    subplot(3,1,[1 2]);
    iLS = rem(idm,lenLS)+1;
%     for j=1:room.lmCount
%         iLC = rem(j,lenLC)+1;
%         iMK = rem(j,lenMK)+1;
%         plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
%         plot(rxX,SNRdb(:,j),plStyle);
%     end
%     iMK = rem(j+1,lenMK)+1;
    iLC = rem(idm,lenLC)+1;
    iMK = rem(idm,lenMK)+1;
    plStyle = [clLC{iLC} clLS{iLS} clMK{iMK}];
    if(fPLOT(idm))
        LGDH(end+1) = plot(rxX,max(SNRdb,[],2),plStyle);
        LGD{end+1} = sprintf('m=%-2.0d, %-2.2f W',txm,room.luminaire.rdFlux);
        %     axis([min(rxX) max(rxX) 25 2*ceil(max(SNRdb(:))/2)]);
        axis([min(rxX) max(rxX) 25 50]);
        grid on;
        % plot illuminance,irradiance
        figure(figSNR);
        subplot(3,3,6+iSF);
        room.drawIlluminance(cLocation(ilX,ilY,ilZ),ctILLTh,gca);
        figure(figSNR);
        subplot(3,3,6+iSF);
        caxis([0 2500]);
        axis([0 room.L 0 room.W 0 2500]);
        title(sprintf('Illumination (m=%d), P_{tx}= %-2.2f W',txm,room.luminaire.rdFlux));
        
        room.drawIlluminance(cLocation(ilX,ilY,ilZ));
        figILL(end+1) = gcf;
        caxis([0 2500]);
        axis([0 room.L 0 room.W 0 2500]);
        title(sprintf('Illumination (m=%d), P_{tx}= %-2.2f W',txm,room.luminaire.rdFlux));
        iSF = iSF + 1;
    end
    ITER = ITER + 1;
end
figure(figSNR);
legend(LGDH,LGD);
if fSAVEALL
    fname = [ctDirRes num2str(ifig) ' SNR' CHARIDXARCHIVE];
    saveas(gcf,[fname '.fig'],'fig');
    saveas(gcf,[fname '.png'],'png');
    ifig = ifig + 1;
    
    for f=1:numel(figILL)
        fname = [ctDirRes num2str(ifig) ' illumination' CHARIDXARCHIVE];
        saveas(figILL(f),[fname '.fig'],'fig');
        saveas(figILL(f),[fname '.png'],'png');
        ifig = ifig + 1;
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





















