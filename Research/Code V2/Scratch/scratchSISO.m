% scratch
close all;
clearvars;
clc;

daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fFILTBLUE = true;
fNOISESHOTONLY = false;

%% VARIABLE PARAMS
rngIllSet = 50:50:800; % varILLSet
rngDim = ((max(rngIllSet)-rngIllSet)*100)./max(rngIllSet); 
idIll400 = find(rngIllSet == 400);
ctDimRef100 = 400;
% rngIllSet = 400; % varILLSet
% rngTxD = [0.5 0.8 1.0 1.5 2];
% rngLkD = 1:0.2:3;
rngLkD = 0.2:0.2:3;
idLk2 = find(rngLkD == 2);
% rngTxa = [1e-3 1e-2 10e-2 20e-2 50e-2];
% rngRxa = [5e-3 1e-2 5e-2 10e-2 20e-2];

%%  Set variables
% varIllSet = 400; % varILLSet
% varLkD = 2;
varTxa = 1e-2;
varRxa = 5e-3;
varRxZeta = pi/6;
varRxAlpha = 0;
varRxTau = 0;
%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\3. s-MIMO\5. CapacitySVDKnownCSI\SVD-VLC full angle\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scratchSISO.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/3. s-MIMO/5. CapacitySVDKnownCSI/SVD-VLC full angle/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scratchSISO.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\3. s-MIMO\\5. CapacitySVDKnownCSI\\SVD-VLC full angle\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scratchSISO.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'SISO.m'];
ctFileVars = [ctDirRes 'SISOdata.mat'];
ctILLTh = 300;
ctrxB = 50e6;

%% constants
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
% Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,17.1335*Wpsd);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);
% room.setIlluminance(.) will set the output power of the transmitters

Ambpsd = 5.8e-6*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

if fFILTBLUE
    % Bf = getSOG(470,10,1,lambdas);
    % Bf = 0.7*Bf/max(Bf);
    Bf = ones(size(lambdas));
    Bf(lambdas > 300 & lambdas < 500) = 1;
end

%% initialize room
room = cRoom(4,4,4);
txLoc = cLocation(room.L/2,room.W/2,3);
room.luminaire = cLuminaire(Wch,txLoc);
% room.setIlluminance(cLocation(room.L/2,room.W/2,1),ctILLSet);

C_LID = zeros(numel(rngLkD),numel(rngIllSet));
for idIllSet = 1:numel(rngIllSet)
    varIllSet = rngIllSet(idIllSet);
    room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
    
    for idLkD = 1:numel(rngLkD)
        varLkD = rngLkD(idLkD);
        
        [rxX, rxY, rxZ] = getGrid(1,1,1,2,2,2,'Fill');
        % [rxX,rxY,rxZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
        % [rxX rxY rxZ] = getGrid(1,1,1,1,1,1);
        rxX = rxX+room.L/2;
        rxY = rxY+room.W/2;
        rxZ = rxZ + 3-varLkD;
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
%         Ignore shot noise
%         SNR = (Isig.^2)./(In+Ish);
        SNR = (Isig.^2)./In;
        C_LID(idLkD,idIllSet) = log2(1 + SNR);
        
    end
end
openfig([ctDirRes 'fig_CapvsDimvsAp']);
figure(gcf);
hold all;
plot(rngDim,C_LID(idLk2,:),'r-.');
[~,~,~,lStr] = legend;
nStr = numel(lStr);
lStr{nStr+1} = 'SISO';
legend(lStr,'Location','NorthWest');
axis auto;
set(gca, 'xdir','reverse');
if fSAVEALL
    fname = [ctDirRes 'fig_CapvsDimvsAp2.png'];
    saveas(gcf,fname,'png');
end


openfig([ctDirRes 'fig_CapvsLinkvsAp']);
figure(gcf);
hold all;
plot(rngLkD,C_LID(:,idIll400),'r-.');
[~,~,~,lStr] = legend;
nStr = numel(lStr);
lStr{nStr+1} = 'SISO';
legend(lStr,'Location','NorthWest');
axis auto;

if fSAVEALL
    fname = [ctDirRes 'fig_CapvsLinkvsAp2.png'];
    saveas(gcf,fname,'png');
end
% SNR = (Isig.^2)./(In+Ish);
% SNRdb = 10*log10(SNR);

% %% FIGURES
% if fSAVEALL
%     delete([ctDirRes '*.png']);
% end
% ifig = 1;
%
% % draw setup
% figure;
% room.drawSetup(cLocation(room.L/2,room.W/2,rxZ(1)),cOrientation(0,0,0),rcvr.rxFOV);
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' Setup.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% rotate3d on;
%
% % plot LED PSDs
% figure;
% plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax);
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Normalized PSD');
% % title('Normalized PSDs of LEDs');
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' LED PSD.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
%
% % CIE plot
% figure;
% obs = cCIE;
% obs.getCoordinates(Wch.npsd);
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' CIEplot.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
%
% % plot filter responses
% figure;
% plot(rcvr.sensor.filter(1).X,rcvr.sensor.filter(1).Y);
% % title('Filter Transmission');
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Filter Transmission');
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' FilterResponse.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
%
% % plot receiver responses
% figure;
% plot(rcvr.sensor.responsivity(1).X,rcvr.sensor.responsivity(1).Y);
% % title('Sensor Responsivity');
% grid on;
% xlabel('Wavelength (nm)');
% ylabel('Sensor Responsivity');
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' SensorResponsivity.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
%
% % plot illuminance,irradiance
% figure;
% H(1) = gca;
% rotate3d on;
% figure;
% H(2) = gca;
% rotate3d on;
% room.drawIlluminance(rcvr.location,ctILLTh,H);
% if fSAVEALL
%     figure(get(H(1),'Parent'));
%     fname = [ctDirRes num2str(ifig) ' Illuminance 3D.png'];
%     saveas(gcf,fname,'png');
%     view(2);
%     fname = [ctDirRes num2str(ifig) ' Illuminance 2D.png'];
%     saveas(gcf,fname,'png');
%     view(3);
%     ifig = ifig + 1;
%
%     figure(get(H(2),'Parent'));
%     fname = [ctDirRes num2str(ifig) ' Illuminance Coverage.png'];
%     saveas(gcf,fname,'png');
%     ifig = ifig + 1;
% end
% figure;
% room.drawIrradiance;
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' Irradiance 3D.png'];
%     saveas(gcf,fname,'png');
%     view(2);
%     fname = [ctDirRes num2str(ifig) ' Irradiance 2D.png'];
%     saveas(gcf,fname,'png');
%     view(3);
%     ifig = ifig + 1;
% end
% rotate3d on;
%
% % plot SNR
% figure;
% surf(rxX,rxY,reshape(SNRdb,size(rxX)),'FaceColor','interp');
% axis([0 room.L 0 room.W min(SNRdb(:)) max(SNRdb(:))]);
% grid on;
% colorbar;
% view(3);
% tStr = sprintf('SNR\nMax: %0.2f dB',max(SNRdb(:)));
% title(tStr);
% xlabel('X');
% ylabel('Y');
% if fSAVEALL
%     fname = [ctDirRes num2str(ifig) ' SNR Tx ' num2str(txLoc.X,'%0.2f') ' ' num2str(txLoc.Y,'%0.2f') ' ' num2str(txLoc.Z,'%0.2f') ' 3D.png'];
%     saveas(gcf,fname,'png');
%     view(2);
%     fname = [ctDirRes num2str(ifig) ' SNR Tx ' num2str(txLoc.X,'%0.2f') ' ' num2str(txLoc.Y,'%0.2f') ' ' num2str(txLoc.Z,'%0.2f') ' 2D.png'];
%     saveas(gcf,fname,'png');
%     view(3);
%     ifig = ifig + 1;
% end
% rotate3d on;
% other
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest);
    save(ctFileVars);
end

set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);























