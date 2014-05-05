close all
clearvars
clc

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP

LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
obs = cCIE;
% [x,y] = obs.getCoordinates(Wch.npsd);
% obs.getCoordinates(Wch.npsd);
x = 0.3823;
y = 0.3837;
%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\BeamSteering\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrBeamSteering.m';
        ctMatDir = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Matfiles\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/BeamSteering/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrBeamSteering.m';
        ctMatDir = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Matfiles/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\MatlabResults\\BeamSteering\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrBeamSteering.m';
        ctMatDir = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Matfiles\\';
    otherwise
        error('Station not defined');
end

RGBledmat = [ctMatDir 'defLEDrgb01.mat'];
if exist(RGBledmat,'file')
    load(RGBledmat,'RGB');
else
    RES = 0.1;
    RGB = cLEDrgb(RES);
    RGB.initialize();
    save(RGBledmat,'RGB');
end

[S,R,G,B,tr,tg,tb] = RGB.getPSD(x,y);
scl = Wch.lmFlux/S.lmFlux;
S.scaleOutputFlux(scl);
R.scaleOutputFlux(scl);
G.scaleOutputFlux(scl);
B.scaleOutputFlux(scl);
CCT = mccamy(x,y);
[x1,y1] = planckXY(CCT);
figure;
obs.getCoordinates(S.npsd);

figure;
subplot(2,1,1);
plot(S.npsd.X,S.npsd.Y);
axis tight;
title(sprintf('white psd\nCCT=%0.1fK',CCT));
subplot(2,3,6);
plot(R.npsd.X,R.npsd.Y);
axis tight;
title('red psd');
subplot(2,3,5);
plot(G.npsd.X,G.npsd.Y);
axis tight;
title('green psd');
subplot(2,3,4);
plot(B.npsd.X,B.npsd.Y);
axis tight;
title('blue psd');