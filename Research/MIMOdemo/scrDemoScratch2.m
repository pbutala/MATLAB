% scrDemoScratch
%*************************************************************************%
%****************************** TRANSMIT *********************************%
%*************************************************************************%
close all;
clearvars;
clc;

% Initialize system parameters
S = SYSTEMTX();
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
txFr = signal(:);
frIP = updnClock(txFr,S.dFs,S.dCLKs,true);

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
ITER = 1;
% Initialize system parameters
S = SYSTEMRX();
Stx = SYSTEMTX();
READSHIFT = 4;
READSKIP = floor(Stx.dCLKs/S.dCLKs); 
READMAX = length(TXFRAME)+READSHIFT-1;
READIdx = rem(READSHIFT:READSKIP:READMAX,length(TXFRAME))+1;
rxFr = TXFRAME(READIdx);
%*************************************************************************%
%****************************** RECEIVE **********************************%
%*************************************************************************%
frmIP = updnClock(rxFr,S.dCLKs,S.dFs,true);
RXFRMSHFT = frmIP(1:length(txFr));

pSTART = floor(READSHIFT*S.dCLKs/Stx.dCLKs)+1;
            
% FRAME RECONSTRUCT
RXFRAME = [RXFRMSHFT(pSTART:end);RXFRMSHFT(1:pSTART-1)];
RXDATA = RXFRAME;
% CALCULATE OFFSET AND GAIN
CHNLOFST = min(RXDATA);
% CHNLGAIN = (max(RXPILOT)-min(RXPILOT))/Stx.dSIGMAX;
CHNLGAIN = (max(RXDATA)-min(RXDATA))/((max(txFr)-min(txFr))*Stx.dSIGMAX);

% SCALE RECEIVED SIGNAL FOR DETECTION
RXFRAMESCL = (RXFRAME-CHNLOFST)/(CHNLGAIN*Stx.dSIGMAX);
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
























