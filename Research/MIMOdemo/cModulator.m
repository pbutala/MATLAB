classdef cModulator
    % base class to handle modulation schemes.
    properties(SetAccess = immutable, GetAccess = protected)
        BUFIN;                      % Input buffer (tx bits)
        BUFOUT;                     % Output buffer (tx signal)
    end % properties - immutable, protected
    
    properties(SetAccess = immutable)
        CLKIN;                      % Clock-In, Data sample clock
        CLKOUT;                     % Clock-Out, Signal sample clock
        SMPLF;                      % Sampling factor (In to Out)
        SIGLO;                      % Modulated signal lo clip
        SIGHI;                      % Modulated signal hi clip
        SCALE;                      % Scale signal (after clip)
        OFFSET;                     % Signal offset (after clip and scale)
    end % properties - immutable
    
    properties(Abstract, SetAccess = protected)
        BPSYM;                      % Bits Per Symbol (input)
        NPSYM;                      % Samples Per Symbol (output)
    end % properties - abstract, protected
    
    methods
        % CONSTRUCTOR
        function obj = cModulator(clkin, clkout, SZIN, SZOUT, lo, hi, scl, ofst)
        % cModulator class constructor
            obj.BUFIN = cFIFO(SZIN);
            if~exist('SZOUT','var')
                SZOUT = ceil(SZIN*clkout/clkin);
            end
            if~exist('lo','var')
                lo = -power(2,15);
            end
            if~exist('hi','var')
                hi = power(2,15)-1;
            end
            if~exist('scl','var')
                scl = 1;
            end
            if~exist('ofst','var')
                ofst = 0;
            end
            obj.BUFOUT = cFIFO(SZOUT);
            obj.CLKIN = clkin;
            obj.CLKOUT = clkout;
            obj.SMPLF = clkout/clkin;
            obj.SIGLO = lo;
            obj.SIGHI = hi;
            obj.SCALE = scl;
            obj.OFFSET = ofst;
        end % contructor
        
        % STREAM INPUT
        function write(obj,data)
            obj.BUFIN.enQ(data);
            sig = obj.modulate();
            sigUS = updnClock(sig,obj.CLKIN,obj.CLKOUT,false);
            sigUS(sigUS<obj.SIGLO) = obj.SIGLO;
            sigUS(sigUS>obj.SIGHI) = obj.SIGHI;
            sigUS = sigUS*obj.SCALE + obj.OFFSET;
            obj.BUFOUT.enQ(sigUS);
        end % write
        
        % STREAM OUTOUT
        function sig = read(obj,MAXOUT)
            sig = obj.BUFOUT.deQ(MAXOUT);
        end % read
    end % methods
    
    methods(Abstract)
        sig = modulate(obj);  % convert bits to signal and queue in BUFOUT.
    end % methods - abstract
end