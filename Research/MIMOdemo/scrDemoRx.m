% Create figure to plot received signal and processing chain
FIGRX = figure('Name', 'Demo Receive', 'NumberTitle', FIGTITLE);
FIGNR = 2;
FIGNC = 2;

% start tcp server on interface
cmd = sprintf('start C:\\ProgramData\\_4DSP_Training\\FMC116\\Debug\\Fmc116APP.exe 1 ML605 %d %d',...
    demo.ADC.dETHID,demo.ADC.dCLKSRC);
dos(cmd);
pause(2);                                                                   % wait till the interface app initializes

% create tcp client
RXCLT=tcpip('localhost', demo.ADC.SKT_PORT, 'NetworkRole', 'client');
% set BurstSize for data to receive
BS = 32*ceil(demo.frmRxNSmp16/32);                                           % BurstSize for FMC116 must be an integer multiple of 32

% Connect to FMC116 interface
RXCLT.InputBufferSize = BS*2;                                               % FMC116 BurstSize*2 (for int16)
fopen(RXCLT);

% set BurstSize for data to receive
fwrite(RXCLT,[demo.ADC.CMD_BURSTSIZE demo.ADC.CHNL_1 demo.ADC.LEN_BS_LSB demo.ADC.LEN_BS_MSB]);
fwrite(RXCLT,typecast(int16(BS),'uint8'));  

% Start acquisition routine
CHNLREAD = demo.ADC.CHNL_1;
fwrite(RXCLT,[demo.ADC.CMD_DATA CHNLREAD demo.ADC.ZERO_UC demo.ADC.ZERO_UC]);      % Read data from channel 1

BYTECOUNT = 0;
while (BYTECOUNT < BS*2)
    BYTECOUNT = RXCLT.BytesAvailable;
end

% READ DATA FROM RECEIVER
frameLong = typecast(uint8(fread(RXCLT,BS*2,'uchar')),'int16');
frame = frameLong(1:demo.frmRxNSmp16);
frame = double(frame*-1);

% Plot receive Frame (@ADC clock)
figure(FIGRX);
subplot(FIGNR,FIGNC,1);
plot(1:demo.frmRxNSmp16, frame);
axis([1 demo.frmRxNSmp16 min(frame) max(frame)]);
xlabel('Normalized time');
ylabel('Signal value');
title('Receive frame (@ADC clock)');

% % Start acquiring next channel
% switch(CHNLREAD)
%     case demo.ADC.CHNL_1
%         CHNLREAD = demo.ADC.CHNL_2;
%         % CHNLREAD = demo.ADC.CHNL_1;
%         SPIDX = 1;
%     case demo.ADC.CHNL_2
%         CHNLREAD = demo.ADC.CHNL_3;
%         % CHNLREAD = demo.ADC.CHNL_1;
%         SPIDX = 2;
%     case demo.ADC.CHNL_3
%         CHNLREAD = demo.ADC.CHNL_4;
%         SPIDX = 3;
%     case demo.ADC.CHNL_4
%         CHNLREAD = demo.ADC.CHNL_1;
%         SPIDX = 4;
% end
% fwrite(RXCLT,[demo.ADC.CMD_DATA CHNLREAD demo.ADC.ZERO_UC demo.ADC.ZERO_UC]);

% Upsample received frame to transmit (DAC) clock for better alignment
% rxFrmUP = updnClock(frame,demo.ADC.dCLKs,demo.DAC.dCLKs);

% Plot receive Frame (@DAC clock)
% figure(FIGRX);
% subplot(FIGNR,FIGNC,2);
% plot(1:demo.frmTxNSmp16, rxFrmUP);
% axis([1 demo.frmTxNSmp16 min(rxFrmUP) max(rxFrmUP)]);
% xlabel('Normalized time');
% ylabel('Signal value');
% title('Receive frame (@DAC clock)');

% Find frame starting index by aligning pilot
% idx = demo.plt.alignPilot(rxFrmUP-min(rxFrmUP), demo.DAC.dCLKs);
idx = demo.plt.alignPilot(frame-min(frame), demo.ADC.dCLKs);

% Find pilot in up-sampled received frame
% NFRM = demo.frmTxNSmp16;
% pltIs = idx:idx+demo.plt.LENGTH-1;
% while ~isempty(pltIs(pltIs<=0))
%     pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
% end
% while ~isempty(pltIs(pltIs>NFRM))
%     pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
% end
% rxPltUS = rxFrmUP(pltIs);

% Find pilot in received frame
NFRM = demo.frmRxNSmp16;
pltIs = idx:idx+demo.plt.getLength(demo.ADC.dCLKs)-1;
while ~isempty(pltIs(pltIs<=0))
    pltIs(pltIs<=0) = pltIs(pltIs<=0) + NFRM;
end
while ~isempty(pltIs(pltIs>NFRM))
    pltIs(pltIs>NFRM) = pltIs(pltIs>NFRM) - NFRM;
end
rxPltUS = frame(pltIs);

% Show pilot start and stop cursors on figure
figure(FIGRX);
subplot(FIGNR,FIGNC,1);
drCursor([pltIs(1), pltIs(end)],'Vertical',['g-';'r-']);

% Estimate channel gain from pilot
% hhat = demo.plt.getScale(rxPltUS, demo.DAC.dCLKs);
hhat = demo.plt.getScale(rxPltUS, demo.ADC.dCLKs);
scl = 1/(hhat*demo.mod.SCALE*demo.DAC.dGAIN);

% Find signal in frame
if pltIs(end) < pltIs(1) % if index roll over
    sigIs = pltIs(end)+1:pltIs(1)-1;
else
    sigIs = [pltIs(end)+1:NFRM, 1:pltIs(1)-1];
end

% Show signal start and stop cursors on figure
figure(FIGRX);
subplot(FIGNR,FIGNC,1);
drCursor([sigIs(1), sigIs(end)],'Vertical',['g:';'r:']);

% rxSigUS = rxFrmUP(sigIs);
% rxSig = updnClock(rxSigUS,demo.DAC.dCLKs,demo.ADC.dCLKs,'IDEALRECT');
% rxSig = updnClock(rxSigUS,demo.DAC.dCLKs,demo.ADC.dCLKs,'RAISEDCOSINE');

rxSig = frame(sigIs);

% Plot receive signal (@fSig*2 clock)
figure(FIGRX);
subplot(FIGNR,FIGNC,3);
plot(1:demo.demod.NPSYM, rxSig);
axis([1 demo.demod.NPSYM min(rxSig) max(rxSig)]);
xlabel('Normalized time');
ylabel('Signal value');
title('Recovered signal');

% Scale received signal
rxSig = rxSig*scl;

% Demodulate received signal
rxSig = rxSig - min(rxSig);
demo.demod.write(rxSig);
rxBits = demo.demod.read(demo.demod.COUNTOUT);

% Find bit errors
berr = biterr2(txBits,rxBits);
fprintf('BERRs = %3d\n',berr);

fclose(RXCLT);
delete(RXCLT);
clear RXCLT;