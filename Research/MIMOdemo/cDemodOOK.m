classdef cDemodOOK < cDemodulator
    % class to handle OOK demodulation
    % 1:ON, 0:OFF
    properties(Constant)
        ON_BIT = 1;                 % Bit level for ON
        OFF_BIT = 0;                % Bit level for OFF
    end % properties - constant
    
    properties(SetAccess = immutable)
        ON_SIG;                     % Signal level for ON
        OFF_SIG;                    % Signal level for OFF
        THRESHOLD;                  % Signal demodulation threshold (ON + OFF)/2
    end % properties - immutable
    
    properties(SetAccess = protected)
        BPSYM;                      % Bits Per Symbol
        NPSYM;                      % Samples Per Symbol
    end % properties - protected
    
    properties
        FILTER = 'IDEALRECT';       % Up/Dn sampling filter type
    end % properties
    
    methods
        % CONSTRUCTOR
        function obj = cDemodOOK(on, off, clkin, clkout, bufszin, bufszout)
            bps = 1;
            nps = ceil(clkin/clkout);
            if~exist('clkin','var')
                clkin = 1;                  % default input clock
            end
            if~exist('clkout','var')
                clkout = 1;                 % default output clock
            end
            if~exist('bufszin','var')
                bufszin = nps;            % default input buffer size
            end
            if~exist('bufszout','var')
                bufszout = bufszin/nps;  % default output buffer size
            end
            obj = obj@cDemodulator(clkin, clkout, bufszin, bufszout);
            obj.ON_SIG = on;
            obj.OFF_SIG = off;
            obj.THRESHOLD = (on + off)/2;
            obj.BPSYM = bps;
            obj.NPSYM = nps;
        end % constructor
    end % methods
    
    methods 
        % DEMODULATE
        function demodulate(obj)
            while(obj.BUFIN.COUNT >= obj.NPSYM)
                strm = obj.BUFIN.deQ(obj.NPSYM);
                % soft decoding (max vote)
                I = (strm > obj.THRESHOLD);
                bit = mode(double(I));
                obj.BUFOUT.enQ(bit);
            end % while loop
        end % modulate
    end % methods - overloaded
end