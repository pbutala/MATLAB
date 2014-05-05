% scrOFDMWDM
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
set(0,'DefaultFigurePaperPosition',[0 0 8 6]);
dlinems = get(0,'DefaultLineMarkerSize');
set(0,'DefaultLineMarkerSize',4);

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fARCHIVE = false;
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = '1_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

% constants
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;
RMN = 627; RSD = 10; RSC = 1;
GMN = 530; GSD = 10; GSC = 1;
BMN = 470; BSD = 10; BSC = 1;
RES = 0.1;
sPSDTYP = 'Gaussian';
TOTALBITS = 1e4;
BERTH = 1e-3;

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\12 WDMOFDM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrOFDMWDM.m';
        ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\LEDPSD\';
        sPSDDIR = [ctMatDir sPSDTYP '\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                   RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\'];
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/12 WDMOFDM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrOFDMWDM.m';
        ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/LEDPSD/';
        sPSDDIR = [ctMatDir sPSDTYP '/' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d',...
                   RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '/'];
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\LEDPSD\\';
        sPSDDIR = [ctMatDir sPSDTYP '\\' sprintf('R_%d_%d_%d_G_%d_%d_%d_B_%d_%d_%d\\',...
                   RMN,RSD,RSC,GMN,GSD,GSC,BMN,BSD,BSC) '\\'];
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrOFDMWDM' CHARIDXARCHIVE '.m'];
ctFileVars = [ctDirRes STRPREFIX 'datOFDMWDM' CHARIDXARCHIVE '.mat'];      % file to store workspace
RGBledmat = [sPSDDIR sprintf('res_%0.5f',RES) '.mat'];

%% config
lkIl = 400;
lkPl = 1;       % plane for illumination
lkTx = 3;       % plane for transmitters
rmL = 4; rmW = 4; rmH = 4; % room L,W,H
rxX = 2; rxY = 2; rxZ = lkPl; % receiver location
flT = 0.7;

%% PSDs

Rpsd = getSOG(RMN,RSD,RSC,lambdas);
Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);

Gpsd = getSOG(GMN,GSD,GSC,lambdas);
Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);

Bpsd = getSOG(BMN,BSD,BSC,lambdas);
Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);

Ambpsd = 5.8e-2*ones(size(lambdas));        % Ambient PSD
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);   % Ambient Channel

%% ranges
rngCCT = 25000;
lenCCT = numel(rngCCT);

%% variables
% TODO: NEED MORE RESOLUTION ON XYZ
if exist(RGBledmat,'file')
    load(RGBledmat,'RGB');
else
    RGB = cLEDrgb(RES,Rch,Gch,Bch);
    RGB.initialize();
    if ~exist(sPSDDIR,'dir')
        mkdir(sPSDDIR);
    end
    save(RGBledmat,'RGB');
end

% initialize receiver filter models
% TODO: May want to use 'Lorentzian' model
Rf = getSOG(RMN,RSD,RSC,lambdas);
Rf = flT*Rf/max(Rf);
Gf = getSOG(GMN,GSD,GSC,lambdas);
Gf = flT*Gf/max(Gf);
Bf = getSOG(BMN,BSD,BSC,lambdas);
Bf = flT*Bf/max(Bf);

% initialize receivers
Rrx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Rrx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rf);
Grx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Grx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gf);
Brx = cSinglePixelReceiverWhiteReflection(rxX,rxY,rxZ);
Brx.sensor.filter = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);

%% logic
BITSTREAM = randi([0 1],[1,TOTALBITS]);
for iT = 1:lenCCT
    % initialize transmitter
    [x,y] = planckXY(rngCCT(iT));
    [~,R,G,B] = RGB.getPSD(x,y);
    
    % set room with new transmitter
    clear room;
    room = cRoom(rmL,rmW,rmH);        % create room and set size (4m x 4m x 4m)
    room.luminaire = cLuminaire([R G B],rmL/2,rmW/2,lkTx);
    room.luminaire.orientation = cOrientation(pi,0,0);
    locCntr = cLocation(rmL/2,rmW/2,lkPl);
    room.setIlluminance(locCntr,lkIl);        % set illuminance at receiver location
    % S = R+G+B; S = room.luminaire.color
    % R,G,B are the luminaire channels
    % lmFluxClr and rdFluxClr are the individual channel fluxes
    
    % calculate H red
    tHfs = room.getFreeSpaceGain(Rrx.location,Rrx.orientation,Rrx.rxFOV);
    [HRr,~] = Rrx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HRg,~] = Rrx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HRb,~] = Rrx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H green
    tHfs = room.getFreeSpaceGain(Grx.location,Grx.orientation,Grx.rxFOV);
    [HGr,~] = Grx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HGg,~] = Grx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HGb,~] = Grx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % calculate H blue
    tHfs = room.getFreeSpaceGain(Brx.location,Brx.orientation,Brx.rxFOV);
    [HBr,~] = Brx.getSignal(R./R.rdFlux,tHfs(1),Ambch);
    [HBg,~] = Brx.getSignal(G./G.rdFlux,tHfs(2),Ambch);
    [HBb,~] = Brx.getSignal(B./B.rdFlux,tHfs(3),Ambch);
    % construct channel matrix (does not include transmitted power)
    H = [HRr HRg HRb;...
         HGr HGg HGb;...
         HBr HBg HBb];
     
    % Monte-Carlo O-OFDM runs
    
end













% restore defaults
set(0,'DefaultLineMarkerSize',dlinems);
set(0,'DefaultLineLineWidth',dlinelw);
set(0,'DefaultAxesFontName',daxesfontname);
% set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);
set(0,'DefaultFigurePaperPosition',dfigpp);
set(0,'DefaultFigurePaperUnits',dfigpu);
set(0,'DefaultFigurePaperPositionMode',dfigppm);

disp(' ');
disp('--DONE--');





































