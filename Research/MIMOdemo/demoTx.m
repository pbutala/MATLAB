% DEMO TRANSMIT

% START FMC204APP AS A LOCAL SERVER
% ENABLE DEMO.M AS CLIENT
% CONNECT
% GENERATE OFDM DATA
% TRANSFER A BLOCK OF DATA TO FMC204APP
% START TRANSMIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clearvars;
clc;

% GLOBAL SETUP
global S TXCLT TXRUN FIGTITLE pilot;

TXRUN = 1;
% Initialize system parameters
S = SYSTEMTX();
FIGTITLE = 'Off';
% Generate Pilot
pilot = S.getPILOT01();

rng('default');
% start tcp server on interface
cmd = sprintf('start C:\\ProgramData\\_4DSP_Training\\FMC204\\Debug\\Fmc204APP.exe 1 ML605 %d %d',S.dETHID,S.dCLKSRC);
dos(cmd);
pause(2); % wait till the interface app initializes

% start tcp client and Connect to FMC204 interface
TXCLT=tcpip('localhost', 30001, 'NetworkRole', 'client');
TXCLT.OutputBufferSize = S.frmLEN8SF;
fopen(TXCLT);
% set BurstSize for data to transmit
fwrite(TXCLT,typecast(int16(S.frmLEN16SF),'uint8'));

% Start trasmission routine
txTmr = timer('Name','StartTxClt','StartDelay', 0.5, 'Period', 1, 'TasksToExecute', 1, ...
    'ExecutionMode', 'fixedRate');
txTmr.StartFcn = {@fnTxTimer};
txTmr.TimerFcn = {@fnTxTimer};
txTmr.StopFcn = {@fnTxTimer};
start(txTmr); 

























