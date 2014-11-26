classdef cPilot
    % base class to handle pilots.
    
    properties(SetAccess = protected)
        PILOT;                      % Pilot samples at CLKOUT rate
    end % properties - protected
    
    properties(SetAccess = immutable, GetAccess = protected)
        CLKIN;                      % Clock-In, Data sample clock
        CLKOUT;                     % Clock-Out, Signal sample clock
    end % properties - immutable, protected
    
    properties(Dependent = true, SetAccess = private)
        LENGTH;                     % Length of output pilot signal
    end % properties - dependent
    
    methods
        % CONSTRUCTOR
        function obj = cPilot(clkin, clkout)
        % cPilot class constructor
        obj.CLKIN = clkin;
        obj.CLKOUT = clkout;
        end % contructor
    end % methods
    
    methods(Abstract, Access = protected)
        val = getPilot(obj, clkout);  % gets pilot signal and stores in PILOT
    end % methods - abstract, protected
    
    methods(Abstract)
        idx = alignPilot(obj, sig, clksmp);  % aligns signal containing pilot sampled at clksmp
    end % methods - abstract
    
    methods
        function scl = getScale(obj, plt, clksmp)
            pltR = obj.getPilot(clksmp);
            pltR = pltR - min(pltR);
            pltR = pltR/max(pltR);
            
            if ~isequal(numel(plt),numel(pltR))
                error('Input pilot and reference have different lengths.');
            end
            plt = plt-min(plt);
            pltR = pltR-min(pltR);
            
            scl = sum(abs(plt))/sum(abs(pltR));
        end % getScale
        
        function len = getLength(obj, clksmp)
            len = ceil(obj.LENGTH*clksmp/obj.CLKOUT);
        end % getLength
    end % methods - public
    
    % Getters/Setters
    methods
        function val = get.LENGTH(obj)
            val = numel(obj.PILOT);
        end % LENGTH
    end % methods - getters/setters
    
end