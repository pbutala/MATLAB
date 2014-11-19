classdef cDemodulator
    % base class to handle modulation schemes.
    properties(SetAccess = immutable, GetAccess = protected)
        BUFIN;                      % Input buffer (tx bits)
        BUFOUT;                     % Output buffer (tx signal)
    end % properties - immutable, protected
    
    properties(SetAccess = immutable)
        CLKIN;                      % Clock-In, Data sample clock
        CLKOUT;                     % Clock-Out, Signal sample clock
        SMPLF;                      % Sampling factor (In to Out)
    end % properties - immutable
    
    properties(Abstract, SetAccess = protected)
        BPSYM;                      % Bits Per Symbol  (input)
        NPSYM;                      % Samples Per Symbol  (output)
    end % properties - abstract, protected
    
    methods
        % CONSTRUCTOR
        function obj = cDemodulator(clkin, clkout, SZIN, SZOUT)
        % cModulator class constructor
            obj.BUFIN = cFIFO(SZIN);
            if~exist('SZOUT','var')
                SZOUT = ceil(SZIN*clkout/clkin);
            end
            obj.BUFOUT = cFIFO(SZOUT);
            obj.CLKIN = clkin;
            obj.CLKOUT = clkout;
            obj.SMPLF = clkout/clkin;
        end % contructor
        
        % STREAM INPUT
        function write(obj,data)
            obj.BUFIN.enQ(data);
            obj.demodulate();
        end % write
        
        % STREAM OUTOUT
        function sig = read(obj,MAXOUT)
            sig = obj.BUFOUT.deQ(MAXOUT);
        end % read
    end % methods
    
    methods(Abstract)
        demodulate(obj);  % convert signal to bits and queue in BUFOUT.
    end % methods - abstract
end