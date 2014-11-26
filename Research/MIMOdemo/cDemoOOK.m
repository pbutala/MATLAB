classdef cDemoOOK < cDemoConfig
    % class to handle demo with OOK modulation
    properties(SetAccess = protected)
        demod;                  % demodulator
        mod;                    % modulator
        plt;                    % pilot
    end % properties - protected
    
    methods
        % CONSTRUCTOR
        function obj = cDemoOOK(Fsig, Fplt, BPFrm)
            % constructor - creates an instance of cDemoOOK
            obj = obj@cDemoConfig(Fsig, Fplt);
            obj.plt = cPilotBarker('BARKER13', obj.fPlt, obj.DAC.dCLKs);
            SIGSCALE = (obj.DAC.dSIGMAX - obj.DAC.dSIGMIN);
            obj.mod = cModOOK(obj.DAC.dSIGMAX, obj.DAC.dSIGMIN, obj.fSig, obj.DAC.dCLKs,...
                SIGSCALE, obj.DAC.dSIGMIN, BPFrm);
            obj.mod.FILTER = 'RAISEDCOSINE';
            obj.demod = cDemodOOK(1, 0, obj.ADC.dCLKs, obj.fSig, BPFrm*ceil(obj.ADC.dCLKs/obj.fSig));
        end % cDemoOOK
    end % methods - public
end