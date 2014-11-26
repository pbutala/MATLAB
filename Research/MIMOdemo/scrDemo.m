% scrDemo
close all;
clearvars;
clc;
rng('Default');

FIGTITLE = 'Off';

fSig = 25e6;
fPlt = 25e6;

% % ----OOK----
% fprintf('--OOK--\n');
% BPFrm = 128;
% demo = cDemoOOK(dFs, dFp, BPFrm);
% %------------

% ----OFDM----
fprintf('--OFDM--\n');
demo = cDemoOFDM(fSig, fPlt);
BPFrm = demo.mod.BPSYM;
%------------
% Generate bits to transmit
txBits = randi([0 1],BPFrm,1);

% Start Transmit Routine
scrDemoTx;

% Start Receive Routine
scrDemoRx;

%--END--