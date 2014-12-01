classdef cDemoOFDM < cDemoConfig
    % class to handle demo with OOK modulation
    properties(SetAccess = protected)
        demod;                  % demodulator
        mod;                    % modulator
        plt;                    % pilot
    end % properties - protected
    
    properties(SetAccess = immutable)
        ofdmTyp;                % type of O-OFDM
        ofst;                   % offset applied to time domain signal
        nsc;                    % # total subcarriers
        msc;                    % M-QAM
        cliph;                  % clip high
        clipl;                  % clip low
        syms;                   % QAM symbols
    end % properties - immutable
    
    methods
        % CONSTRUCTOR
        function obj = cDemoOFDM(Fsig, Fplt, spFrm)
            % constructor - creates an instance of cDemoOFDM
            if ~exist('spFrm','var')
                spFrm = 1;
            end
            obj = obj@cDemoConfig(Fsig, Fplt, spFrm);
            obj.plt = cPilotBarker('BARKER13', obj.fPlt, obj.DAC.dCLKs);
            
            % OFDM parameters
            obj.ofdmTyp = 'DCOOFDM';
            obj.ofst = 3.2;
            obj.nsc = 8;
            obj.msc = 4;
            obj.cliph = 3.2;
            obj.clipl = 0;
            obj.syms = getQAMsyms(obj.msc,true);
            
            %        %******************* TO SET SIGNAL SCALE *******************%
            %        [cdf,bins,lo,hi] = getOFDMdist(ofdmTyp, nsc, syms, ofst, 0);
            %        I = find(cdf>0.99,1,'first');
            %        SIGMAX = bins(I);
                     %***********************************************************%
                    % SIGMAX | OFDM | NSC  | MSC  | cdf>TH |
                    % 5.5102 |  DCO |  64  |  16  |   0.99 |
                    %        |  DCO |   8  |   4  |   0.99 |
            SIGMAX = 10;
            SIGMIN = 0;
            SIGSCALE = (obj.DAC.dSIGMAX - obj.DAC.dSIGMIN)/(SIGMAX-SIGMIN);
            obj.mod = cModOFDM(obj.ofdmTyp, obj.nsc, obj.syms, obj.ofst, SIGMIN, SIGMAX,...
                SIGSCALE,...
                obj.DAC.dSIGMIN - SIGMIN,...
                obj.fSig, obj.DAC.dCLKs);
            obj.mod.FILTER = 'IDEALRECT';
            obj.demod = cDemodOFDM(obj.ofdmTyp, obj.nsc, obj.syms, obj.ofst, obj.ADC.dCLKs, obj.fSig);
        end % cDemoOOK
    end % methods - public
end