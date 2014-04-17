% imaging OSM.
clc
close all
clear all

daxesfontname = get(0,'DefaultAxesFontName');
set(0,'DefaultAxesFontName','Helvetica');
daxesfontsize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',16);
dfigvis = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','On');

%%%% Constants %%%%
LAMBDAMIN = 200;
LAMBDADELTA = 1;
LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

s1=18; m1=450; a1=1; s2=60; m2=555; a2=2.15*a1; s3=25; m3=483; a3=-0.2*a1;
Wpsd = getSOG([m1 m2 m3],[s1 s2 s3],[a1 a2 a3],lambdas);
Wpsd = Wpsd/(sum(Wpsd)*LAMBDADELTA);
Wch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Wpsd);

Ambpsd = 5.8e-2*ones(size(lambdas));
Ambch = cPSD(LAMBDAMIN,LAMBDADELTA,LAMBDAMAX,Ambpsd);

%%%% OSM %%%%
m=2; % bits per symbols
M = 2^m; % # of symbols
alphabit = modem.pammod('M', 2^m); %PAM modulation

Nt=4; % # 0f Txs 
Nr = 4; % # of Rxs
snr =0:5:30; % signal-to-noise ratio
loop =1;

%%%% VARIABLE PARAMS
% rngIllSet = 50:50:800; % varILLSet
rngIllSet = 400; % varILLSet
rngDim = ((max(rngIllSet)-rngIllSet)*100)./max(rngIllSet);
idIll400 = find(rngIllSet == 400);
ctDimRef100 = 400;
% rngIllSet = 400; % varILLSet
% rngTxD = [0.5 0.8 1.0 1.5 2];
% rngLkD = 1:0.2:3;
% rngLkD = 0.2:0.2:3;
% rngLkD = [0.9 1.2 1.9 2.2];
rngLkD = 2;
idLk2 = find(rngLkD == 2);
% rngTxa = [1e-3 1e-2 10e-2 20e-2 50e-2];
% rngOpf = [1e-3 5e-3 1e-2 5e-2 10e-2];
% rngOpD = [1e-3 5e-3 1e-2];
% rngOpD = [1e-3 2e-3 (5e-3*2/sqrt(pi))];
rngOpD = 5e-3;
% rngOpD = [1e-3 5e-3];
% rngRxa = [5e-3 1e-2 5e-2 10e-2 20e-2];
% rngRxal = [0.5e-3 1e-3 5e-3 1e-2 5e-2 10e-2];% depends on rngRxa only use values smaller than 'varRxa'
rotD = 15;
rngRxZeta = 0:rotD*pi/180:pi/2-rotD*pi/180;
% rngRxAlpha = -pi:15*pi/180:pi;
rngRxAlpha = 0:rotD*pi/180:pi-rotD*pi/180;

%%%%  Set variables
% varIllSet = 400; % varILLSet
varLkD = 2;
varTxa = 5e-2;
varOpf = 5e-3;
varOpD = 5e-3;
varOpFOV = pi/3;
varRxa = 2e-3;
varRxal = 1e-3;% depends on rngRxa only use values smaller than 'varRxa'
varTxD = 0.2; %D between TX
% varTxD = sqrt(2)*varRxal*(varLkD-varOpf)/varOpf;
varRxZeta = 0;
varRxAlpha = 0;
varRxTau = 0;

fN = varOpf/varOpD;

%%%% Room %%%%
room = cRoom(4,4,4);

% [txX, txY, txZ] = getGrid(room.L,room.W,1,varTxD,varTxD,2,'Fill');
[txX, txY, txZ] = getGrid(2,2,1,varTxD,varTxD,1);
txX = txX + room.L/2;
txY = txY + room.W/2;
txZ = txZ + 3;
txLoc = cLocation(txX,txY,txZ);
txOri = cOrientation(pi,0,0);
room.luminaire = cLuminaire(Wch,txLoc);
room.luminaire.orientation = txOri;
room.luminaire.dimension = cSize(varTxa,varTxa,1e-3);

[plX, plY, plZ] = getGrid(room.L,room.W,1,0.2,0.2,2,'Fill');
plX = plX + room.L/2;
plY = plY + room.W/2;
plZ = plZ + 1;
plLoc = cLocation(plX,plY,plZ);

NPxX = varRxa/varRxal;
NPxY = varRxa/varRxal;
NPxD = varRxal;
% NPxD = 5.91e-2/NPxX;
[pxX, pxY, pxZ] = getGrid(NPxX,NPxY,1,NPxD,NPxD,2);

% [rxX,rxY,rxZ] = getGrid(room.L,room.W,1,1,1,2,'Fill');
[rxX,rxY,rxZ] = getGrid(1,1,1,1,1,1);
rxX = rxX + room.L/2;
rxY = rxY + room.W/2;
rxZ = rxZ + 3 - varLkD;
rxLoc = cLocation(rxX,rxY,rxZ);
rxOri = cOrientation(varRxZeta,varRxAlpha,varRxTau);
Wrx = cImagingReceiverWhiteReflection(rxLoc.X,rxLoc.Y,rxLoc.Z,pxX,pxY,varOpf,fN);

Wrx.orientation = rxOri;
Wrx.optics.fN = fN;
Wrx.optics.T = 1; % set lens transmission to 100%
Wrx.sensor.dimension.L(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.W(1:NPxX,1:NPxY) = NPxD;
Wrx.sensor.dimension.H(1:NPxX,1:NPxY) = 1e-3;
Wrx.optics.fov = varOpFOV;
Wrx.sensor.color = [1 1 1];

%%%%%%%%Imaging Rx Channel Martix%%%%%%%%%%%%%%%%
H = room.getFreeSpaceGain(Wrx.location,Wrx.orientation,Wrx.rxFOV);
nL = numel(room.luminaire);
txSz = room.lmDimension;
% currently assuming just 1 lum group at multiple loc
color = room.luminaire(1).color;
txOri = room.luminaire(1).orientation;
Hs = reshape(H,Wrx.rxCount,sum(room.lmCount(:)));
[Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp, Hch] = ...
    Wrx.getSignal(color,Hs,Ambch,txLoc,txSz,txOri);
% [Ipx,Iambpx,Isig,Iamb,frSpOnPx,frPxOnSp] = ...
%     Wrx.getSignal(color.*(color.lmFlux),Hs,Ambch,txLoc,txSz,txOri);
Hch1 = squeeze(Hch(1,:,:))';
Hch1t = Hch1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj =snr
    err = 0; % Detected Errors
    iter = 1000; % # of Iterations
    for ii = 1:iter
        % bit_stream = randint(1, log2(M*Nt)); % generate bit stream e.g. 4 bits per symbol for 4 Txs
        bit_stream = randi([0,1], [1,log2(M*Nt)]);
        sym_index = bin2dec(int2str(bit_stream(1:log2(M)))); % Symbol index
        led_index = bin2dec(int2str(bit_stream(log2(M)+1:log2(Nt*M))))+1; % LED index
        tx = zeros(Nt,1);
        tx(led_index)=modulate(alphabit,sym_index); % assign Symbol to LED 
        
        %%%%%%%Add AWGN Channel%%%%%%%%%%%%%%%%%%%%%%%%%%
        y= awgn((Hch1t*tx),jj,'measured');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%Optimal Spatial Modulation Detector%%%%%%%%%%%%%%
        for i=1:Nt
            for j=1:M
                g_jq = Hch1t(:,i)*modulate(alphabit,j-1);
                min_2(i,j) = norm(g_jq, 'fro')^2 - 2*real(y'*g_jq); %This vector contains the errors to all antennas.
            end
        end
        [tmp,idx] =  min( min_2(:)) ;
        [idx1,idx2] = find(min_2==tmp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%Estimated bits from LED index%%%%%%%%%%%%%%%%%%%
        tmp1 = dec2bin(idx1-1,log2(Nt));
        
        for uu = 1:log2(Nt)
            led_bits(uu)= bin2dec(tmp1(uu));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%Estimated bits from Symbol%%%%%%%%%%%%%%%%%%%%%%
        tmp2 = dec2bin(idx2-1,log2(M));
        
        for xx = 1:log2(M)
            sym_bits(xx)= bin2dec(tmp2(xx));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rx_bits = [sym_bits led_bits];

        err = err+biterr(rx_bits,bit_stream);
    end
    
        bit_err(loop) = err/((log2(M*Nt))*iter)
        loop = loop+1
end

save OSM_IMRx_4QAM.mat bit_err snr;
figure
hold on
semilogy(snr,bit_err,'-bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5)
xlabel('SNR (dB)')
ylabel('Bit Error Ratio')
title('Imaging Rx OSM 4x4')
grid
hold off

set(0,'DefaultAxesFontName',daxesfontname);
set(0,'DefaultAxesFontSize',daxesfontsize);
set(0,'DefaultFigureVisible',dfigvis);