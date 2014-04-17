% DEMO RECEIVE

% START FMC116APP AS LOCAL SERVER
% ENABLE DEMORX.M AS A CLIENT
% CONNECT
% START RECEIVE
% TRANSFER A BLOCK OF DATA FROM FMC204APP
% RECOVER OFDM DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;

% GLOBAL SETUP
global S Stx RXCLT RXRUN RXFRAME RXDATASCL;
RXRUN = 1;
RXDATASCL = [];
RXFRAME = []; % will be initialized on first read
% Initialize system parameters
S = SYSTEMRX();
Stx = SYSTEMTX();
% Generate Data that was transmitted
rng('default');

% start tcp server on interface
cmd = sprintf('start C:\\ProgramData\\_4DSP_Training\\FMC116\\Debug\\Fmc116APP.exe 1 ML605 %d %d',S.dETHID,S.dCLKSRC);
dos(cmd);
pause(2); % wait till the interface app initializes

% create tcp client
RXCLT=tcpip('localhost', 30002, 'NetworkRole', 'client');
% set BurstSize for data to receive
% BurstSize for FMC116 must be an integer multiple of 32
BS = 32*ceil(S.frmLEN16SF/32);
% Connect to FMC116 interface
RXCLT.InputBufferSize = BS*2; % FMC116 BurstSize*2 (for int16)
fopen(RXCLT);
fwrite(RXCLT,typecast(int16(BS),'uint8'));  

% Start acquisition routine
rxTmr = timer('Name','StartRxClt','StartDelay', 0.5, 'Period', 0.25, 'TasksToExecute', 80, ...
    'ExecutionMode', 'fixedRate');
rxTmr.StartFcn = {@fnRxTimer};
rxTmr.TimerFcn = {@fnRxTimer};
rxTmr.StopFcn = {@fnRxTimer};
start(rxTmr); 































