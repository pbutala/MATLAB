close all;
clearvars;
clc;

% GLOBAL SETUP
S = SYSTEMTX();
USF = S.dSF;

% get signal
sig = S.getPILOT01();
% take fourier transform of signal
sL = S.frmPILOTLEN16;
sNFFT = 2^nextpow2(sL);
sFT = fft(sig,sNFFT)/sL;
sf = S.dFs/2*linspace(0,1,sNFFT/2+1);
sFTM = abs(sFT);

% upsample
sigUS = upsample(sig,USF);
% take fourier transform of upsampled signal
uL = length(sigUS);
uNFFT = 2^nextpow2(uL);
uFT = USF*fft(sigUS,uNFFT)/uL;
uf = S.dCLKs/2*linspace(0,1,uNFFT/2+1);
uFTM = abs(uFT);

% create filter
flt = zeros(uNFFT,1);
flt(uf<=S.dFs/2) = 1;
flt(uf>S.dCLKs-S.dFs/2) = 1;

% filter upsampled signal
iFT = uFT.*flt;
iFTM = abs(iFT);
iL = uL;
iNIFFT = 2^nextpow2(iL);
sigIP = iL*ifft(iFT,iNIFFT,1,'symmetric');
txSig = sigIP(1:USF*sL);

% DRAW FIGS
FIGSIG = figure('Name', 'Signal');
subplot(2,3,1);
plot(sig);
axis tight;

subplot(2,3,4);
plot(sf*1e-6,2*sFTM(1:sNFFT/2+1));
axis tight;
title('Single sided amplitude spectrum of transmit signal');
xlabel('Frequency (MHz)');
ylabel('|sFT|');

subplot(2,3,2);
plot(sigUS);
axis tight;

subplot(2,3,5);
plotyy(uf*1e-6,2*uFTM(1:uNFFT/2+1),uf*1e-6,flt(1:uNFFT/2+1));
axis tight;
title('Single sided amplitude spectrum of upsampled signal');
xlabel('Frequency (MHz)');
ylabel('|uFT|');

subplot(2,3,3);
% plotyy(1:iNIFFT,sigIP,1:USF*sL,txSig);
plot(txSig);
axis tight;

subplot(2,3,6);
plot(uf*1e-6,2*iFTM(1:iNIFFT/2+1));
title('Single sided amplitude spectrum of filtered upsampled signal');
xlabel('Frequency (MHz)');
ylabel('|iFT|');
axis tight;


































