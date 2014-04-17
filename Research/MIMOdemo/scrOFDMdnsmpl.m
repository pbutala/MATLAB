close all;
clearvars;
clc;

% GLOBAL SETUP
Stx = SYSTEMTX();
% get signal
sigTx = Stx.getPILOT01();
frmTx = updnclock(sigTx,Stx.dFs,Stx.dCLKs,true);
% figure;
% plot(frmTx);

Srx = SYSTEMRX();
frmRx = frmTx(round(1:Stx.dCLKs/Srx.dCLKs:end));
sigRx = updnclock(frmRx,Srx.dCLKs,Srx.dFs,true);

% % take fourier transform of signal
% sL = Stx.frmPILOTLEN16;
% sNFFT = 2^nextpow2(sL);
% sFT = fft(sig,sNFFT)/sL;
% sf = Stx.dFs/2*linspace(0,1,sNFFT/2+1);
% sFTM = abs(sFT);
% 
% % upsample
% sigUS = upsample(sig,USF);
% % take fourier transform of upsampled signal
% uL = length(sigUS);
% uNFFT = 2^nextpow2(uL);
% uFT = USF*fft(sigUS,uNFFT)/uL;
% uf = Stx.dCLKs/2*linspace(0,1,uNFFT/2+1);
% uFTM = abs(uFT);
% 
% % create filter
% flt = zeros(uNFFT,1);
% flt(uf<=Stx.dFs/2) = 1;
% flt(uf>Stx.dCLKs-Stx.dFs/2) = 1;
% 
% % filter upsampled signal
% iFT = uFT.*flt;
% iFTM = abs(iFT);
% iL = uL;
% iNIFFT = 2^nextpow2(iL);
% sigIP = iL*ifft(iFT,iNIFFT,1,'symmetric');
% txSig = sigIP(1:USF*sL);
% 
% % DRAW FIGS
% FIGSIG = figure('Name', 'Transmit');
% subplot(2,3,1);
% plot(sig);
% axis tight;
% 
% subplot(2,3,4);
% plot(sf*1e-6,2*sFTM(1:sNFFT/2+1));
% axis tight;
% title('Single sided amplitude spectrum of transmit signal');
% xlabel('Frequency (MHz)');
% ylabel('|sFT|');
% 
% subplot(2,3,2);
% plot(sigUS);
% axis tight;
% 
% subplot(2,3,5);
% ax = plotyy(uf*1e-6,2*uFTM(1:uNFFT/2+1),uf*1e-6,flt(1:uNFFT/2+1));
% axis(ax, 'tight');
% title('Single sided amplitude spectrum of upsampled signal');
% xlabel('Frequency (MHz)');
% ylabel('|uFT|');
% 
% subplot(2,3,3);
% plot(txSig);
% axis tight;
% 
% subplot(2,3,6);
% plot(uf*1e-6,2*iFTM(1:iNIFFT/2+1));
% title('Single sided amplitude spectrum of filtered upsampled signal');
% xlabel('Frequency (MHz)');
% ylabel('|iFT|');
% axis tight;
% 
% % CHANNEL
% clearvars -except txSig;
% % tx clock to rx closk is factor of 8.
% sig = txSig(1:8:end);
% clear txSig;
% % DOWNSAMPLE
% Srx = SYSTEMRX();
% [USF,DSF] = rat(1/Srx.dSF);
% 
% % take fourier transform of signal
% sL = Srx.frmPILOTLEN16SF;
% sNFFT = 2^nextpow2(sL);
% sFT = fft(sig,sNFFT)/sL;
% sf = Srx.dCLKs/2*linspace(0,1,sNFFT/2+1);
% sFTM = abs(sFT);
% 
% % upsample
% sigUS = upsample(sig,USF);
% % take fourier transform of upsampled signal
% uL = length(sigUS);
% uNFFT = 2^nextpow2(uL);
% uFT = USF*fft(sigUS,uNFFT)/uL;
% uf = (Srx.dCLKs*USF)/2*linspace(0,1,uNFFT/2+1);
% uFTM = abs(uFT);
% 
% % create filter
% flt = zeros(uNFFT,1);
% flt(uf<=Srx.dFs/2) = 1;
% flt(uf>Srx.dCLKs*USF-Srx.dFs/2) = 1;
% 
% % filter upsampled signal
% iFT = uFT.*flt;
% iFTM = abs(iFT);
% iL = uL;
% iNIFFT = 2^nextpow2(iL);
% sigIP = iL*ifft(iFT,iNIFFT,1,'symmetric');
% rxSigUS = sigIP(1:USF*sL);
% 
% % downsample signal
% rxSig = rxSigUS(1:DSF:end);
% % take fourier transform of dnsampled signal
% dL = length(rxSig);
% dNFFT = 2^nextpow2(dL);
% dFT = fft(rxSig,dNFFT)/dL;
% df = (Srx.dFs)/2*linspace(0,1,dNFFT/2+1);
% dFTM = abs(dFT);
% 
% % DRAW FIGS
% FIGSIG = figure('Name', 'Receive');
% subplot(2,4,1);
% plot(sig);
% axis tight;
% 
% subplot(2,4,5);
% plot(sf*1e-6,2*sFTM(1:sNFFT/2+1));
% axis tight;
% title('Single sided amplitude spectrum of receive signal');
% xlabel('Frequency (MHz)');
% ylabel('|sFT|');
% 
% subplot(2,4,2);
% plot(sigUS);
% axis tight;
% 
% subplot(2,4,6);
% ax = plotyy(uf*1e-6,2*uFTM(1:uNFFT/2+1),uf*1e-6,flt(1:uNFFT/2+1));
% axis(ax, 'tight');
% title('Single sided amplitude spectrum of upsampled signal');
% xlabel('Frequency (MHz)');
% ylabel('|uFT|');
% 
% subplot(2,4,3);
% plot(rxSigUS);
% axis tight;
% 
% subplot(2,4,7);
% plot(uf*1e-6,2*iFTM(1:iNIFFT/2+1));
% title('Single sided amplitude spectrum of filtered upsampled signal');
% xlabel('Frequency (MHz)');
% ylabel('|iFT|');
% axis tight;
% 
% subplot(2,4,4);
% plot(rxSig);
% axis tight;
% 
% subplot(2,4,8);
% plot(df*1e-6,2*dFTM(1:dNFFT/2+1));
% axis tight;
% title('Single sided amplitude spectrum of downsampled signal');
% xlabel('Frequency (MHz)');
% ylabel('|dFT|');


















