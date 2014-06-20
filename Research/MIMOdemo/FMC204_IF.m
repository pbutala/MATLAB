classdef FMC204_IF
    properties (Constant)
        % Socket
        SKT_PORT = 30001;
        ZERO_UC = uint8(hex2dec('00'));
        
        % Commands
        CMD_BURSTSIZE = uint8(hex2dec('10'));
        CMD_DATA = uint8(hex2dec('20'));
        CMD_ENCHNL = uint8(hex2dec('30'));
        CMD_ARMDAC = uint8(hex2dec('40'));
        
        % ADC Channel
        CHNL_1 = uint8(hex2dec('01'));
        CHNL_2 = uint8(hex2dec('02'));
        CHNL_3 = uint8(hex2dec('04'));
        CHNL_4 = uint8(hex2dec('08'));
        CHNL_ALL = uint8(hex2dec('0F'));
        
        % Lengths
        LEN_BS_MSB = uint8(hex2dec('00'));
        LEN_BS_LSB = uint8(hex2dec('02'));
    end
end