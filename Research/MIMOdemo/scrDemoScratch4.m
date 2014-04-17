% scrDemoScratch
%*************************************************************************%
%****************************** TRANSMIT *********************************%
%*************************************************************************%
close all;
clearvars;
clc;

% Initialize system parameters
S = SYSTEMTX();
% Generate Pilot
pilot = S.getPILOT01();

% Generate Data
rng('default');
data = randi(S.modM,[S.frmDATALEN16 1]);
datab = dec2binMat(data-1,log2(S.modM));

% Generate OFDM symbol
txSig = genOFDMsignal(... % Variable Arguments to the function
    'data',data,...
    'OFDMtype',S.modTYPE,...
    'N',S.modN,...
    'Symbols',getQAMsyms(S.modM,true),...
    'ClipLow',S.modCLIPL,...
    'ClipHigh',S.modCLIPH,...
    'OffsetDcoStddev', S.modOFST,...
    'ShowConst',false);

% Scale OFDM symbol
signal = txSig/S.modCLIPH;  % signal max is 1

% Generate frame
txFr = [pilot(:);signal(:)];
frIP = updnClock(txFr,S.dFs,S.dCLKs,true);

% Ensure Pilot is between 0 and 1
frIP = frIP/(max(frIP(1:S.frmPILOTLEN16SF))- min(frIP(1:S.frmPILOTLEN16SF)));
frIP = frIP - min(frIP(1:S.frmPILOTLEN16SF));

% Shift, Scale
TXFRAME = frIP*S.dSIGMAX;
TXFRAME(TXFRAME > S.dSIGMAX) = S.dSIGMAX;
TXFRAME(TXFRAME < S.dSIGMIN) = S.dSIGMIN;

figure;
ax(1) = subplot(2,1,1);
ax(2) = subplot(2,1,2);
plotFreq(txFr,S.dFs,'single',ax);
axes(ax(1));
hold on;
axes(ax(2));
hold on;

%*************************************************************************%
%*************************** IDEAL CHANNEL********************************%
%*************************************************************************%
clearvars -except TXFRAME txFr data datab ax;
% Initialize system parameters
S = SYSTEMRX();
Stx = SYSTEMTX();
ITER = 1;

READSHIFT = 15;
READSKIP = floor(Stx.dCLKs/S.dCLKs); 
READMAX = length(TXFRAME)+READSHIFT-1;
READIdx = rem(READSHIFT:READSKIP:READMAX,length(TXFRAME))+1;
rxFr = TXFRAME(READIdx);
%*************************************************************************%
%****************************** RECEIVE **********************************%
%*************************************************************************%

%**************************************
rxFrUP = updnClock(rxFr,S.dCLKs,Stx.dCLKs,false);
RXFRMSHFT = rxFrUP(1:Stx.frmLEN16SF);
frmcyc = [RXFRMSHFT; RXFRMSHFT(1:floor(Stx.frmPILOTLEN16SF))];
frmcyc = frmcyc - min(frmcyc);
Pilot = updnClock(S.getPILOT01(),S.dFs,Stx.dCLKs,false);
[acor, lag] = xcorr(Pilot,frmcyc);
[~,I] = max(acor(lag<=0));
pSTART = rem(abs(lag(I)),Stx.frmLEN16SF)+1;
% FRAME RECONSTRUCT
rxFrUP0 = [RXFRMSHFT(pSTART:end);RXFRMSHFT(1:pSTART-1)];
%**************************************
RXFRAME = updnClock(rxFrUP0,Stx.dCLKs,S.dFs,false);
RXPILOT = RXFRAME(1:S.frmPILOTLEN16);
RXDATA = (RXFRAME(S.frmPILOTLEN16+1:end));
% plotFreq(RXFRMSHFT,S.dCLKs);
% append

            

% CALCULATE OFFSET AND GAIN
CHNLOFST = min(RXPILOT);
% CHNLGAIN = (max(RXPILOT)-min(RXPILOT))/Stx.dSIGMAX;
CHNLGAIN = (max(RXPILOT)-min(RXPILOT))/((max(Pilot)-min(Pilot))*Stx.dSIGMAX);

% SCALE RECEIVED SIGNAL FOR DETECTION
RXFRAMESCL = (RXFRAME-CHNLOFST)/(CHNLGAIN*Stx.dSIGMAX);
% RXDATASCL = S.modCLIPL + ((RXDATA-CHNLOFST)*(S.modCLIPH-S.modCLIPL))/(CHNLGAIN*Stx.dSIGMAX);
RXDATASCL = S.modCLIPL + ((RXDATA-CHNLOFST)*(S.modCLIPH-S.modCLIPL))/(CHNLGAIN*Stx.dSIGMAX);
            
plotFreq(RXFRAMESCL,S.dFs,'Single',ax,'r:');
% DECODE OFDM
rxDat = decodeOFDMsignal(RXDATASCL,...
    'OFDMtype',S.modTYPE,...
    'N',S.modN,...
    'Symbols',getQAMsyms(S.modM,true),...
    'ShowRcv',true);

rxDatb = dec2binMat(rxDat-1,log2(S.modM));
% CALCULATE BIT-ERRORS
err = biterr2(datab,rxDatb);
ber = err/numel(rxDatb);
fprintf('%d. BER = %0.2e\n',ITER,ber);

figure;
subplot(4,1,1);
plot(rxFr);
axis tight;
title('Received frame');
subplot(4,1,2);
plot(rxFrUP);
axis tight;
title('Received frame upsampled');

subplot(4,1,3);
plot(lag,acor);
axis tight;
title(sprintf('Autocorrelation with pilot\nStart Idx: %d',pSTART));
subplot(4,1,4);
plot(RXFRAME);
axis tight;
title('Frame un-shifted');

figure(3);





















