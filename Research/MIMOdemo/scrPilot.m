% pilot filtering
close all;
clearvars;
clc;

% Initialize system parameters
S = SYSTEMTX();

% Generate Pilot
P = S.getPILOT01();
[Fb,Fa] = butter(4,0.8);
Pf = filter(Fb,Fa,P);
freqz(Fb,Fa);

% pilots frequency response
L = length(P);
N = 2^nextpow2(L);
FT = fft(P,N)/L;
f = (S.dFs/2)*linspace(0,1,N/2+1);
FTM = abs(FT);

% filtered pilots frequency response
fFT = fft(Pf,N)/L;
fFTM = abs(fFT);

% FIGURES
figure;
subplot(2,2,1);
plot(P);

subplot(2,2,3);
plot(f*1e-6,2*FTM(1:N/2+1));

subplot(2,2,2);
plot(Pf);

subplot(2,2,4);
plot(f*1e-6,2*fFTM(1:N/2+1));