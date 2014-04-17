% scratch
close all;
clearvars;
clc;

%% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = true;
fFILTBLUE = false;
fNOISESHOTONLY = false;

%% VARIABLE PARAMS
rngIllSet = 200:50:800; % varILLSet
% rngIllSet = 400; % varILLSet
% rngTxD = [0.5 0.8 1.0 1.5 2];
% rngLkD = 1:1:3;
% rngTxa = [1e-3 1e-2 10e-2 20e-2 50e-2];
% rngOpf = [1e-3 5e-3 1e-2 5e-2 10e-2];
rngOpD = [1e-3 5e-3 1e-2 5e-2 10e-2 20e-2];
% rngRxa = [5e-3 1e-2 5e-2 10e-2 20e-2];
% rngRxal = [0.5e-3 1e-3 5e-3 1e-2 5e-2 10e-2];% depends on rngRxa only use values smaller than 'varRxa'

%%  Set variables
% varIllSet = 400; % varILLSet
varLkD = 2;
varTxa = 1e-2;
varOpf = 5e-3;
% varOpD = 5e-3;
varOpFOV = pi/3;
varRxa = 18e-3;
varRxal = 1e-3;% depends on rngRxa only use values smaller than 'varRxa'
varTxD = sqrt(2)*varRxal*(varLkD-varOpf)/varOpf;

mtSNRvLxvApMAX = zeros(numel(rngIllSet),numel(rngOpD));

%% SETUP
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\3. s-MIMO\4. TxSpacing\2. PointSource Iamp\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrTxSpacing.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/3. s-MIMO/4. TxSpacing/2. PointSource Iamp/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrTxSpacing.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\3. s-MIMO\\4. TxSpacing\\2. PointSource Iamp\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrTxSpacing.m';
    otherwise
        error('Station not defined');
end

ctFileCodeDest = [ctDirRes 'scrTxSpacing.m'];
ctIllTh = 300;

%% constants
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


%% initialize room
room = cRoom(7,7,4);
% [txX txY txZ] = getGrid(2,2,1,varTxD,varTxD,2);
[txX, txY, txZ] = getGrid(room.L,room.W,1,varTxD,varTxD,2,'Fill');
txX = txX + room.L/2;
txY = txY + room.W/2;
txZ = txZ + 3;
txLoc = cLocation(txX,txY,txZ);
txOri = cOrientation(pi,0,0);
room.luminaire = cLuminaire(Wch,txLoc);
room.luminaire.orientation = txOri;
room.luminaire.dimension = cSize(varTxa,varTxa,1e-3);

[rxX, rxY, rxZ] = getGrid(1,1,1,2,2,2,'Fill');
% [rxX rxY rxZ] = getGrid(room.L,room.W,1,1,1,2,'Fill');
% [rxX rxY rxZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
rxX = rxX + room.L/2;
rxY = rxY + room.W/2;
rxZ = rxZ + 3 - varLkD;
rxLoc = cLocation(rxX,rxY,rxZ);

[plX, plY, plZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
plX = plX + room.L/2;
plY = plY + room.W/2;
plZ = plZ + 3 - varLkD;
plLoc = cLocation(plX,plY,plZ);

NPxX = varRxa/varRxal;
NPxY = varRxa/varRxal;
NPxD = varRxal;
% NPxD = 5.91e-2/NPxX;
[pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);
Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,varOpf,1);
Wrx.optics.T = 1; % set lens transmission to 100%
Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;
Wrx.optics.fov = varOpFOV;

if fFILTBLUE
    Wrx.sensor.filter(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Bf);
    Wrx.sensor.responsivity(1:NPxX,1:NPxY) = cCurve(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,0.1*ones(1,LAMBDAMAX-LAMBDAMIN+1));
end

for idxOpD = 1:numel(rngOpD)
    varOpD = rngOpD(idxOpD);
    for idxIllSet = 1:numel(rngIllSet)
        varIllSet = rngIllSet(idxIllSet);
        fN = varOpf/varOpD;
        Wrx.optics.fN = fN;
        
        flPreFix = sprintf('%0.0flx_TxD_%0.1fcm_LkD_%0.0fm_Txa_%0.0fmm_Opf_%0.0fmm_OpD_%0.0fmm_Rxa_%0.0fmm_NpxXY_%0.0f__',...
            varIllSet,varTxD*1e2,varLkD,varTxa*1e3,varOpf*1e3,varOpD*1e3,varRxa*1e3,varRxa/varRxal);
        ctFileVars = [ctDirRes flPreFix 'txSpacing.mat'];
        %% Logic
        disp(flPreFix);
        
        room.setIlluminance(cLocation(room.L/2,room.W/2,1),varIllSet);
        
        % for each receiver group
        disp('Calculating Channel Matrix');
        H = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);
        
        nL = numel(room.luminaire);
        txSz = room.lmDimension;
        % currently assuming just 1 lum group at multiple loc
        color = room.luminaire(1).color;
        txOri = room.luminaire(1).orientation;
        Hs = reshape(H,Wrx.rxCount,sum(room.lmCount(:)));
        disp('Calculating Signal and Noise');
        [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = ...
            Wrx.getSignal(color,Hs,Ambch,txLoc,txSz,txOri);
        
        if fNOISESHOTONLY
            In = 0;
        else
            In = Wrx.tia.getNoise(40e6);
        end
        Ish = 2*1.6e-19*(Isig+Iamb)*40e6;
        SNR = (Isig.^2)./(In(1)+Ish);
        SNRdb = 10*log10(SNR);
        
        mtSNRvLxvApMAX(idxIllSet,idxOpD) = max(SNRdb(:));
        
        cLim = [floor(min(SNRdb(:))) ceil(max(SNRdb(:)))];
        
        disp('Saving Results');
        %% SAVE workspace
        if fSAVEALL
            %                                     copyfile(ctFileCodeSrc,ctFileCodeDest);
            save(ctFileVars);
        end
        
        locCntr = cLocation(room.L/2,room.W/2,rxZ(1));
        
        %% FIGURES
        if fSAVEALL
            %                                     delete([ctDirRes '*.png']);
        end
        ifig = 1;
        
        % draw spots
        % loc2 = cLocation(room.L/2+0.25,room.W/2,rxZ(1));
        % figure('Visible','Off');
        figure;
        h1 = gca;
        % h2 = subplot(1,2,2);
        % Wrx.getImage(locCntr+loc2,txLoc,txSz,txOri,[h1 h2]);
        Wrx.getImage(locCntr,txLoc,txSz,txOri,h1);
        if fSAVEALL
            fname = [ctDirRes flPreFix num2str(ifig) ' Spots.png'];
            saveas(gcf,fname,'png');
            ifig = ifig + 1;
        end
        if fCLOSEFIGS
            close;
        end
        % plot illuminance,irradiance
        % [plX plY plZ] = getGrid(room.L,room.W,1,0.5,0.5,2,'Fill');
        % plLoc = cLocation(plX+room.L/2,plY+room.W/2,plZ+0.85);
        % figure('Visible','Off');
        figure;
        H(1) = gca;
        rotate3d on;
        figure;
        H(2) = gca;
        rotate3d on;
        room.drawIlluminance(plLoc,ctIllTh,H);
        if fSAVEALL
            figure(get(H(1),'Parent'));
            fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 3D.png'];
            saveas(gcf,fname,'png');
            view(2);
            fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance 2D.png'];
            saveas(gcf,fname,'png');
            view(3);
            if fCLOSEFIGS
                close;
            end
            ifig = ifig + 1;
            figure(get(H(2),'Parent'));
            fname = [ctDirRes flPreFix num2str(ifig) ' Illuminance Coverage.png'];
            saveas(gcf,fname,'png');
            ifig = ifig + 1;
            if fCLOSEFIGS
                close;
            end
        end
        % figure('Visible','Off');
        figure;
        room.drawIrradiance(plLoc);
        if fSAVEALL
            fname = [ctDirRes flPreFix num2str(ifig) ' Irradiance 3D.png'];
            saveas(gcf,fname,'png');
            view(2);
            fname = [ctDirRes flPreFix num2str(ifig) ' Irradiance 2D.png'];
            saveas(gcf,fname,'png');
            view(3);
            ifig = ifig + 1;
        end
        if fCLOSEFIGS
            close;
        else
            rotate3d on;
        end
        
        %     % plot SNRS
        %     % for each transmitter calculate SNR mesh and plot it
        %     [nRx nTx] = size(SNRdb);
        %     for iT = 1:nTx
        %         snrdb = SNRdb(:,iT);
        %         % figure('Visible','Off');
        %         figure;
        %         set(gca,'climmode','manual');
        %         set(gca,'CLim',cLim);
        %         surf(rxX,rxY,reshape(snrdb,size(rxX)),'FaceColor','interp');
        %         axis([0 room.L 0 room.W min(snrdb(:)) max(snrdb(:))]);
        %         grid on;
        %         colorbar;
        %         view(3);
        %         tStr = sprintf('SNR for Tx @ (%0.2f %0.2f %0.2f)\nMax: %0.2f dB',txLoc.X(iT),txLoc.Y(iT),txLoc.Z(iT),max(snrdb(:)));
        %         title(tStr);
        %         xlabel('X');
        %         ylabel('Y');
        %         if fSAVEALL
        %             fname = [ctDirRes flPreFix num2str(ifig) ' SNR Tx ' num2str(txLoc.X(iT),'%0.2f') ' ' num2str(txLoc.Y(iT),'%0.2f') ' ' num2str(txLoc.Z(iT),'%0.2f') ' 3D.png'];
        %             saveas(gcf,fname,'png');
        %             view(2);
        %             fname = [ctDirRes flPreFix num2str(ifig) ' SNR Tx ' num2str(txLoc.X(iT),'%0.2f') ' ' num2str(txLoc.Y(iT),'%0.2f') ' ' num2str(txLoc.Z(iT),'%0.2f') ' 2D.png'];
        %             saveas(gcf,fname,'png');
        %             view(3);
        %             ifig = ifig + 1;
        %         end
        %         if fCLOSEFIGS
        %             close;
        %         else
        %             rotate3d on;
        %         end
        %     end
        
    end
end
% draw setup
% figure('Visible','Off');
figure;
room.drawSetup(locCntr,cOrientation(0,0,0),Wrx.rxFOV);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' Setup.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
else
    rotate3d on;
end

% plot LED PSDs
% figure('Visible','Off');
figure;
plot(Wch.npsd.X,Wch.npsd.Y/Wch.npsd.Ymax);
grid on;
xlabel('Wavelength (nm)');
ylabel('Normalized PSD');
% title('Normalized PSDs of LEDs');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' LED PSD.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% CIE plot
% figure('Visible','Off');
figure;
obs = cCIE;
obs.getCoordinates(Wch.npsd);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' CIEplot.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% plot filter responses
% figure('Visible','Off');
figure;
plot(Wrx.sensor.filter(1).X,Wrx.sensor.filter(1).Y);
title('Filter response');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' FilterResponse.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end
% plot receiver responses
% figure('Visible','Off');
figure;
plot(Wrx.sensor.responsivity(1).X,Wrx.sensor.responsivity(1).Y);
% title('Sensor Responsivity');
grid on;
xlabel('Wavelength (nm)');
ylabel('Sensor Responsivity');
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' SensorResponsivity.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

% Draw SNR vs Ill vs Aperture plot
figure;
plot(rngIllSet,mtSNRvLxvApMAX(:,1));
lgd{1} = num2str(rngOpD(1));
hold all;
for idxOpD = 2:numel(rngOpD)
    varOpD = rngOpD(idxOpD);
    plot(rngIllSet,mtSNRvLxvApMAX(:,idxOpD));
    lgd{idxOpD} = num2str(rngOpD(idxOpD));
end
grid on;
xlabel('Illumination (lx)');
ylabel('SNR (dB)');
legend(lgd);
if fSAVEALL
    fname = [ctDirRes flPreFix num2str(ifig) ' SNR v Ill v Ap.png'];
    saveas(gcf,fname,'png');
    ifig = ifig + 1;
end
if fCLOSEFIGS
    close;
end

















































