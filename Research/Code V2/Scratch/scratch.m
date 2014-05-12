close all;
clearvars;
clc;

%% FLAGS
fSTATION = 4;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus

LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

RMN = 627; RSD = 10; RSC = 1;               % Mean, SD and scale to generate SPD of Red led
GMN = 530; GSD = 10; GSC = 1;               % Mean, SD and scale to generate SPD of Green led
BMN = 470; BSD = 10; BSC = 1;               % Mean, SD and scale to generate SPD of Blue led
RES = 0.1;                                  % x,y Resolution for xy<->CCT conversion
sSPDTYP = 'Gaussian';                       % LED SPD model

RNGCCT = [3000:3000:6000];                         LENCCT = numel(RNGCCT);              % CCT 
RNGCCTPL = RNGCCT;                     LENCCTPL = numel(RNGCCTPL);          % CCTs to plot

obs = cCIE;
%% SETUP
switch fSTATION
%     case 1
%         ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\BeamSteering\';
%         ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrBeamSteering.m';
%         ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\';
%     case 2
%         ctDirRes = '/home/pbutala/My Documents/MatlabResults/BeamSteering/';
%         ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrBeamSteering.m';
%         ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/';
%     case 3
%         ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\MatlabResults\\BeamSteering\\';
%         ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrBeamSteering.m';
%         ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\';
    case 4
        ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\12 WDMOFDM\\';
        ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrOFDMWDM.m';
        ctMatDir = 'C:\\Users\\Pankil\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\';
    otherwise
        error('Station not defined');
end

%% SPDs
Rpsd = getSOG(RMN,RSD,RSC,lambdas);                 
Rpsd = Rpsd/(sum(Rpsd)*LAMBDADELTA);
Rch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Rpsd);   % Normalized RED SPD

Gpsd = getSOG(GMN,GSD,GSC,lambdas);
Gpsd = Gpsd/(sum(Gpsd)*LAMBDADELTA);
Gch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Gpsd);   % Normalized GREEN SPD

Bpsd = getSOG(BMN,BSD,BSC,lambdas);
Bpsd = Bpsd/(sum(Bpsd)*LAMBDADELTA);
Bch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bpsd);   % Normalized BLUE SPD

RGBledmat = [ctMatDir 'temp~' '.mat'];
if exist(RGBledmat,'file')
    load(RGBledmat,'RGB');
else
    RGB = cLEDrgb(RES,Rch,Gch,Bch);
    RGB.initialize();
    if ~exist(ctMatDir,'dir')
        mkdir(ctMatDir);                 % Create directory to store characterization data
    end
    save(RGBledmat,'RGB');
end

PLACOLC = 'm'; PLACOLS = '-'; PLACOMK = 'o';
PLDCOLC = 'c'; PLDCOLS = '-'; PLDCOMK = 'd';
PLTXLCS = {'r';'g';'b'}; PLTXLSS = {'--';'-.';':'}; PLTXMKS = {'>';'s';'*'};
% Figure CCT config
FIGCCT = figure('Name',sprintf('SPD vs CCT'),'NumberTitle','OFF');
FIGCCTPLNR = power(2,floor(log2(LENCCTPL)/2));
FIGCCTPLNC = power(2,ceil(log2(LENCCTPL)/2));
FIGLDAMIN = 400; FIGLDAMAX = 800;
iTPL = 1;
for iT = 1:LENCCT                                                           % LOOP START CCT
    if RNGCCT(iT) == RNGCCTPL(iTPL)                                         % If CCT selected for plot
        figure(FIGCCT);
        subplot(FIGCCTPLNR,FIGCCTPLNC,iTPL);
        [x,y] = planckXY(RNGCCT(iT)); 
%         fprintf('x=%0.2f y=%0.2f\n',x,y);
        [S,R,G,B,tr,tg,tb] = RGB.getPSD(x,y);                               % Get PSDs at CCT
        plot(R.npsd.X,(tr/S.npsd.Ymax)*R.npsd.Y,PLTXLCS{1});                % Plot RED SPD
        hold on;
        plot(G.npsd.X,(tg/S.npsd.Ymax)*G.npsd.Y,PLTXLCS{2});                % Plot Green SPD
        plot(B.npsd.X,(tb/S.npsd.Ymax)*B.npsd.Y,PLTXLCS{3});                % Plot Blue SPD
        axis([FIGLDAMIN FIGLDAMAX 0 1]);
        xlabel('Wavelength (nm)');
        ylabel('Normalized SPD');
        title(sprintf('CCT = %dK, [x,y] = [%0.2f,%0.2f]',RNGCCTPL(iTPL),x,y));
        iTPL = iTPL + 1;
        figure;
        obs.getCoordinates(S.npsd);
    end    
end