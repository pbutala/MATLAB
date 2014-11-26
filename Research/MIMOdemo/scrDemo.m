% scrDemo
close all;
clearvars;
clc;
rng('Default');

global FIGTITLE demo txBits;
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
% scrDemoTx;
txTmr = timer('Name','StartTxClt',...
              'StartDelay', 4,...
              'Period', 1,...
              'TasksToExecute', 1, ...
              'ExecutionMode', 'fixedRate',...
              'StartFcn', @demoTxTimer,...
              'StopFcn', @demoTxTimer,...
              'TimerFcn', @demoTxTimer);
start(txTmr);

pause(2);

% Start Receive Routine
% scrDemoRx;
rxTmr = timer('Name','StartRxClt',...
              'StartDelay', 4,...
              'Period', 1,...
              'TasksToExecute', 4, ...
              'ExecutionMode', 'fixedRate',...
              'StartFcn', @demoRxTimer,...
              'StopFcn', @demoRxTimer,...
              'TimerFcn', @demoRxTimer);
start(rxTmr);
%--END--