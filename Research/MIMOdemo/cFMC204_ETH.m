classdef cFMC204_ETH
    properties (Constant)
        % Socket
        SKT_PORT = 30001;                       % Ethernet Socket Port
        ZERO_UC = uint8(hex2dec('00'));         % unsigned char ZERO
        
        % Commands
        CMD_BURSTSIZE = uint8(hex2dec('10'));   % Command: Burst sise
        CMD_DATA = uint8(hex2dec('20'));        % Command: Data
        CMD_ENCHNL = uint8(hex2dec('30'));      % Command: Enable Channel(s)
        CMD_ARMDAC = uint8(hex2dec('40'));      % Command: Arm FMC204 DAC
        
        % Channels
        CHNL_1 = uint8(hex2dec('01'));          % Channel: 1
        CHNL_2 = uint8(hex2dec('02'));          % Channel: 2
        CHNL_3 = uint8(hex2dec('04'));          % Channel: 3
        CHNL_4 = uint8(hex2dec('08'));          % Channel: 4
        CHNL_ALL = uint8(hex2dec('0F'));        % Channel: All
        
        % Lengths
        LEN_BS_MSB = uint8(hex2dec('00'));      % Length: BurstSize MSB
        LEN_BS_LSB = uint8(hex2dec('02'));      % Length: BurstSize LSB
        
        % Device specification
        dCLKs = 1e9;                            % device transmit sample clock
        dCHNLCNT = 4;                           % number of channels
        dBITS = 16;                             % Number of bits
        dSIGMAX = power(2,15)-1;                % device max input
        dSIGMIN = -power(2,15);                 % device min input
        dOUTMAX = 0.5;                          % device max output (0.5V)
        dOUTMIN = -0.5;                         % device min output (-0.5V)
        
    end % properties - constant
    
    properties(SetAccess = private)
        dETHID = cFMCETHID.GetNetPCIcards('TX');    % device network card id
        dCLKSRC = 0;                                % device clock source for interface
        dGAIN;                                      % device gain (o/i)
    end % properties private
    
    methods
        % constructor
        function obj = cFMC204_ETH()
            obj.dGAIN = (obj.dOUTMAX - obj.dOUTMIN)/(obj.dSIGMAX - obj.dSIGMIN);   
        end % cFMC204_ETH
        
        % Transfer Function
        function output = getOutput(obj,input)
            output = obj.dGAIN*input;
            output(output > obj.dOUTMAX) = obj.dOUTMAX;
            output(output < obj.dOUTMIN) = obj.dOUTMIN;
        end % getOutput
        
        function val = setRail2Rail(obj,sig)
            lo = min(sig);
            hi = max(sig);
            if(isequal(hi,lo))
                dhi = abs(sig(1) - obj.dSIGMAX);
                dlo = abs(sig(1) - obj.dSIGMIN);
                if(dhi < dlo)
                    val = obj.dSIGMAX*ones(size(sig));
                elseif(dlo < dhi)
                    val = obj.dSIGMIN*ones(size(sig));
                else
                    val = zeros(size(sig));
                end
            else
                slp = (obj.dSIGMAX - obj.dSIGMIN)/(hi-lo);
                val = obj.dSIGMIN + (sig-lo)*slp;
            end
        end % setRail2Rail
        
    end % methods
    
end