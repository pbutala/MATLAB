classdef cDemoConfig
    % class to handle demo configuration
    properties(Abstract = true, SetAccess = protected)
        demod;                  % demodulator
        mod;                    % modulator
        plt;                    % pilot
    end % properties - abstract
    
    properties(SetAccess = immutable)
        ADC;                    % Analog to Digital convertor
        DAC;                    % Digital to Analog convertor
        fPlt;                   % Pilot bandwidth
        fSig;                   % Signal bandwidth
        SPFRM;                  % Number of symbols per frame
    end % properties - immutable
    
    properties(Dependent = true, SetAccess = private)
        frmTxNSmp16;            % number of 16 bit-samples in upsampled frame tx
        frmTxNSmp8;             % number of 8 bit-samples in upsampled frame tx
        frmRxNSmp16;            % number of 16 bit-samples in upsampled frame rx
        frmRxNSmp8;             % number of 8 bit-samples in upsampled frame rx
    end % properties - dependent
    
    methods
        % CONSTRUCTOR
        function obj = cDemoConfig(Fsig, Fplt, spFrm)
            obj.ADC = cFMC116_ETH();
            obj.DAC = cFMC204_ETH();
            obj.fSig = Fsig;
            obj.fPlt = Fplt;
            obj.SPFRM = spFrm;
        end % cDemoConfig
    end % methods
    
    % Getters/Setters
    methods
        function val = get.frmTxNSmp16(obj)
            val = (obj.plt.LENGTH + obj.SPFRM*obj.mod.NPSYM);
        end % frmTxNSmp16
        
        function val = get.frmRxNSmp16(obj)
            val = ceil(obj.frmTxNSmp16*obj.ADC.dCLKs/obj.DAC.dCLKs);
        end % frmRxNSmp16
        
        function val = get.frmTxNSmp8(obj)
            val = obj.frmTxNSmp16*2;
        end % frmTxNSmp8
        
        function val = get.frmRxNSmp8(obj)
            val = obj.frmRxNSmp16*2;
        end % frmRxNSmp8
        
    end % methods (getters/setters)
end