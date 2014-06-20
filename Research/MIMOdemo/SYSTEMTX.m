classdef SYSTEMTX < SYSTEM
    properties(SetAccess = private)
        dCLKs = 1e9;                    % device transmit/receive sample clock
        dCLKSRC = 0;                    % device clock source for interface
        dSIGMAX = power(2,15)-1;        % transmitter peak signal
        dSIGMIN = 0;
        dETHID = SYSTEM.GetNetPCIcards('TX'); % device network card id
    end
end








































