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
            warning('Not yet reliable. Error check needed.');
            obj = obj@cPilot(clkin,clkout);
            obj.NCYCLE = ncycl;
            obj.PILOT = obj.getPilot(clkout);
        end % contructor
    end % methods
    
    methods(Access = protected)
        function val = getPilot(obj, clkout)
            n = 0:ceil(obj.NCYCLE*clkout/obj.CLKIN);
            val = sin(2*pi*n(:)*obj.CLKIN/clkout)-1;
            val = updnClock(val, clkout, clkout, obj.FILTER, false);
        end % getPilot
    end % methods - protected
    
    methods
        function idx = alignPilot(obj, sig, clksmp)  % aligns signal containing pilot sampled at clksmp
            pltR = obj.getPilot(clksmp);
            pltR = pltR - min(pltR);
            pltRLen = numel(pltR);
            sigLen = numel(sig);
            if(sigLen < pltRLen)
                warning('Signal length is smaller than pilot length. Pilot cannot be aligned');
                idx = 1; 
            else
                cycl = sig(1:pltRLen);
                sigC = [sig(:); cycl(:)];
                sigC = sigC - min(sigC);
                [acor, lag] = xcorr(pltR,sigC);
                [~,I] = max(acor(lag<=0));
                idx = rem(abs(lag(I)), sigLen) + 1;
                idx = obj.alignFine(sig, pltR, idx);
            end
        end % alignPilot
    end % methods
    
    methods(Access = private)
        function idx = alignFine(obj,sig, plt, idx0)
            sigLen = numel(sig);
            pltLen = numel(plt);
            DIDX = 3;
            IDXMIN = idx0-DIDX;
            IDXMAX = idx0+DIDX;
            
            MSE = zeros(2*DIDX+1,1);
            IDXs = IDXMIN:IDXMAX;
            IDXs(IDXs<=0) = IDXs(IDXs<=0) + sigLen;
            IDXs(IDXs>sigLen) = IDXs(IDXs>sigLen) - sigLen;
            
            for i=1:numel(IDXs)
                idx = IDXs(i);
                pltI = idx:idx+pltLen-1;
                while ~isempty(pltI(pltI<=0))
                    pltI(pltI<=0) = pltI(pltI<=0) + sigLen;
                end
                while ~isempty(pltI(pltI>sigLen))
                    pltI(pltI>sigLen) = pltI(pltI>sigLen) - sigLen;
                end
                pltR = sig(pltI);
                MSE(i) = sum(plt.*pltR);
            end
            idx = IDXs(MSE==max(MSE));
        end % alignFine
    end % methods - private
end