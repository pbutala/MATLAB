classdef SYSTEMRX < SYSTEM
    properties(SetAccess = private)
        % RECEIVER
        dCLKs = 125e6;                   % Rx 125 Ms/s
%         dCLKs = 1e9;                   % Rx 125 Ms/s
        dCLKSRC = 0;                   % 0: internal, 1: external, 2: internal with ext reference
        dETHID = SYSTEM.GetNetPCIcards('RX'); % device network card id
    end
end








































