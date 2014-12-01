classdef cPilotTone < cPilot
    % base class to handle tone pilots.
    properties(SetAccess = immutable)
        NCYCLE;
    end % properties - immutable
    
    properties
        FILTER = 'IDEALRECT';    % Up/Dn sampling filter type
    end % properties
    
    methods
        % CONSTRUCTOR
        function obj = cPilotTone(ncycl, clkin, clkout)
            % cPilotTone class constructor
            obj = obj@cPilot(clkin,clkout);
            obj.NCYCLE = ncycl;
            obj.PILOT = obj.getPilot(clkout);
        end % contructor
    end % methods
    
    methods(Access = protected)
        function val = getPilot(obj, clkout)
            n = 0:ceil(obj.NCYCLE*clkout/obj.CLKIN);
            val = sin(2*pi*n(:)*obj.CLKIN/clkout);
%             val = updnClock(val, clkout, clkout, obj.FILTER, false);
        end % getPilot
    end % methods - protected
    
    methods
        function idx = alignPilot(obj, sig, clksmp)  % aligns signal containing pilot sampled at clksmp
            pltR = obj.getPilot(clksmp);
            pltR = pltR - min(pltR);
            pltRLen = numel(pltR);
            sigLen = numel(sig);
            MSE = realmax('double')*ones(sigLen,1);
            if(sigLen < pltRLen)
                warning('Signal length is smaller than pilot length. Pilot cannot be aligned');
                idx = 1; 
            else
            scl = 2/(max(sig) - min(sig));
            sig = sig*scl;
            sig = sig - min(sig);
            sig = [sig;sig(1:pltRLen-1)];
            for i=1:sigLen
                seg = sig(i:i+pltRLen-1);
                d = pltR-seg;
                MSE(i) = sum(d.*d);
            end
            idx = find(MSE == min(MSE),1,'first');
            end
        end % alignPilot
    end % methods
end