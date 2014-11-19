classdef cFMC116_ETH
    properties (Constant)
        % Socket
        SKT_PORT = 30002;                       % Ethernet Socket Port
        ZERO_UC = uint8(hex2dec('00'));         % unsigned char ZERO
        
        % Commands
        CMD_BURSTSIZE = uint8(hex2dec('10'));   % Command: Burst sise
        CMD_DATA = uint8(hex2dec('20'));        % Command: Data
        
        % ADC Channel
        CHNL_1 = uint8(hex2dec('01'));          % Channel: 1
        CHNL_2 = uint8(hex2dec('02'));          % Channel: 2
        CHNL_3 = uint8(hex2dec('04'));          % Channel: 3
        CHNL_4 = uint8(hex2dec('08'));          % Channel: 4
        
        % Lengths
        LEN_BS_MSB = uint8(hex2dec('00'));      % Length: BurstSize MSB
        LEN_BS_LSB = uint8(hex2dec('02'));      % Length: BurstSize LSB
        
        % Device specification
        dCLKs = 125e6;                          % device receive sample clock
        dCHNLCNT = 16;                          % number of channels
        dBITS = 14;                             % Channel resolution
        dINP2P= 2;                              % Input P-P V range (2Vpp)
        dOUTMAX = power(2,13)-1;                % device max input
        dOUTMIN = -power(2,13);                 % device min input
        dINMAX = 1;                             % device max input
        dINMIN = -1;                            % device min input
        
    end % properties - constant
    
    properties(SetAccess = private)
        dETHID = cFMCETHID.GetNetPCIcards('RX');    % device network card id
        dCLKSRC = 0;                                % 0: internal, 1: external, 2: internal with ext reference
        dGAIN;                                      % device gain (o/i)
    end % properties private
    
    methods
        % constructor
        function obj = cFMC116_ETH()
            obj.dGAIN = (obj.dOUTMAX - obj.dOUTMIN)/(obj.dINMAX - obj.dINMIN);   
        end % cFMC116_ETH
        
        % Transfer Function
        function output = getOutput(obj,input)
            % assume 0 offset
            output = obj.dGAIN*input;
            output(output > obj.dOUTMAX) = obj.dOUTMAX;
            output(output < obj.dOUTMIN) = obj.dOUTMIN;
        end % getOutput
    end % methods
end