close all;
clearvars;
clc;

% Initialize system parameters
S = SYSTEMTX();

% Generate Data
rng('default');
% data = randi(S.modM,[S.frmDATALEN16 1]);
data = 2*ones(S.frmDATALEN16, 1);

% Generate OFDM symbol
txSig = genOFDMsignal(... % Variable Arguments to the function
    'data',data,...          
    'OFDMtype',S.modTYPE,...    
    'N',S.modN,...
    'Symbols',getQAMsyms(S.modM),...
    'ClipLow',S.modCLIPL,...
    'ClipHigh',S.modCLIPH,...
    'OffsetDcoStddev', S.modOFST);   

% % take fourier transform of signal
% L = length(txSig);
% N = 2^nextpow2(L);
% FT = fft(txSig,N)/L;
% f = S.dFs/2*linspace(0,1,N/2+1);
% FTM = abs(FT);
% FTA = angle(FT);

figure('Name','frequency analysis of OFDM signal','NumberTitle','Off');
subplot(3,1,1);
plot(txSig);
xlabel('Sample');

ax(1) = subplot(3,1,2);
ax(2) = subplot(3,1,3);
plotFreq(txSig,S.dFs,'single',ax);

% subplot(3,1,2);
% plot(f*1e-6,FTM(1:N/2+1));
% xlabel('Frequency (MHz)');
% ylabel('|X_f|');
% title('Single sided amplitude spectrum (DC set to 0)');
% 
% subplot(3,1,3);
% plot(f*1e-6,180*FTA(1:N/2+1)/S.ctPI);
% xlabel('Frequency (MHz)');
% ylabel('|\theta_f^o|');
% title('Single sided angle spectrum');
