classdef cModOOK < cModulator
    % class to handle OOK modulation
    % 1:ON, 0:OFF
    properties(Constant)
        ON_BIT = 1;                 % Bit level for ON
        OFF_BIT = 0;                % Bit level for OFF
    end % properties - constant
    
    properties(SetAccess = protected)
        BPSYM;                      % Bits Per Symbol (input)
        NPSYM;                      % Samples Per Symbol (output)
    end % properties - protected
    
    properties
        FILTER = 'RAISEDCOSINE';    % Up/Dn sampling filter type
    end % properties
    
    methods
        % CONSTRUCTOR
        function obj = cModOOK(on, off, clkin, clkout, scl, ofst,...
                bufszin, bufszout)
            bps = 1;                        % bits per symbol
            nps = ceil(clkout/clkin);       % samples per symbol
            if~exist('clkin','var')
                clkin = 1;                  % default input clock
            end
            if~exist('clkout','var')
                clkout = 1;                 % default output clock
            end
            if~exist('bufszin','var')
                bufszin = bps;            % default input buffer size
            end
            if~exist('bufszout','var')
                bufszout = bufszin*nps;  % default output buffer size
            end
            obj = obj@cModulator(clkin, clkout, bufszin, bufszout, off, on, scl, ofst);
            obj.BPSYM = bps;
            obj.NPSYM = nps;
        end % constructor
    end % methods
    
    methods 
        % MODULATE
        function sig = modulate(obj)
            tSig = obj.BUFIN.deQ(obj.BUFIN.COUNT);
            I = tSig((tSig~=obj.ON_BIT)&(tSig~=obj.OFF_BIT));
            if ~isempty(I)
                warning('Input stream contains values other than ON(%d) and OFF(%d) bits',obj.ON_BIT,obj.OFF_BIT);
            end
            % sig = updnClock(tSig,obj.CLKIN,obj.CLKOUT,obj.FILTER,false);
            sigUS = repmat(tSig(:).',ceil(obj.CLKOUT/obj.CLKIN),1);
            sig = sigUS(:);
        end % modulate
        
    end % methods - overloaded
end