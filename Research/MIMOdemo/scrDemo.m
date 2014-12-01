% scrDemo
close all;
clearvars;
clc;
rng('Default');

global FIGTITLE demo txBits;
FIGTITLE = 'Off';
fSig = 25e6;
fPlt = 25e6;

% ----OOK----
fprintf('--OOK--\n');
spFrm = 256;
demo = cDemoOOK(fSig, fPlt, spFrm);
BPFrm = spFrm;
%------------

% % ----OFDM----
% fprintf('--OFDM--\n');
% spFrm = 1;
% demo = cDemoOFDM(fSig, fPlt, spFrm);
% BPFrm = spFrm*demo.mod.BPSYM;
% %------------

% Generate bits to transmit
txBits = randi([0 1],BPFrm,1);

% Start Transmit Routine
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