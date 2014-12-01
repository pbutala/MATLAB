% scrDemo
close all;
clearvars;
clc;
rng('Default');

global FIGTITLE demo BPFrm datBits;
global SYNC;

FIGTITLE = 'Off';
fSig = 25e6;
fPlt = 25e6;
SYNC = 1;       % 1: transmit 0: receive

% ----OOK----
fprintf('--OOK--\n');
spFrm = 8;
demo = cDemoOOK(fSig, fPlt, spFrm);
BPFrm = spFrm;
%------------

% % ----OFDM----
% fprintf('--OFDM--\n');
% spFrm = 1;
% demo = cDemoOFDM(fSig, fPlt, spFrm);
% BPFrm = spFrm*demo.mod.BPSYM;
% %------------

% initialize buffer to hold data bits for BER calculation
datBits(1) = cFIFO(BPFrm);
datBits(2) = cFIFO(BPFrm);
datBits(3) = cFIFO(BPFrm);
datBits(4) = cFIFO(BPFrm);
% for i=2:4
%     datBits(i) = cFIFO(BPFrm);
% end
%--------------------------------

% Start Transmit Routine
txTmr = timer('Name','StartTxClt',...
              'StartDelay', 2,...
              'Period', 0.25,...
              'TasksToExecute', 7, ...
              'ExecutionMode', 'fixedRate',...
              'StartFcn', @demoTxTimer,...
              'StopFcn', @demoTxTimer,...
              'TimerFcn', @demoTxTimer);
start(txTmr);
%--------------------------------

% Start Receive Routine
rxTmr = timer('Name','StartRxClt',...
              'StartDelay', 2,...
              'Period', 0.25,...
              'TasksToExecute', 7,...
              'ExecutionMode', 'fixedRate',...
              'StartFcn', @demoRxTimer,...
              'StopFcn', @demoRxTimer,...
              'TimerFcn', @demoRxTimer);
start(rxTmr);
%--------------------------------

%--END--