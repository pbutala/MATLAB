classdef cPilotBarker < cPilot
    % base class to handle Barker pilots.
    properties(SetAccess = immutable)
        TYPE;
    end % properties - immutable
    
    methods
        % CONSTRUCTOR
        function obj = cPilotBarker(type, clkin, clkout)
            % cPilotBarker class constructor
            obj = obj@cPilot(clkin,clkout);
            switch(upper(type))
                case {'BARKER11', 'BARKER13'}
                    obj.TYPE = upper(type);
                otherwise
                    error('Pilot type must be ''BARKER11'',''BARKER13''');
            end
            obj.PILOT = obj.getPilot(clkout);
        end % contructor
    end % methods
    
    methods(Access = protected)
        function val = getPilot(obj, clkout)
            switch(obj.TYPE)
                case 'BARKER11'
                    plt = [1;1;1;0;0;0;1;0;0;1;0];
                case 'BARKER13'
                    plt = [1;1;1;1;1;0;0;1;1;0;1;0;1];
            end
            val = updnClock(plt, obj.CLKIN, clkout, true);
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
                idx = 1; scl = 1; ofst = 0;
            else
                cycl = sig(1:pltRLen);
                sigC = [sig(:); cycl(:)];
                sigC = sigC - min(sigC);
                [acor, lag] = xcorr(pltR,sigC);
                [~,I] = max(acor(lag<=0));
                idx = rem(abs(lag(I)), sigLen) + 1;
%                 idxe = idx + pltRLen - 1;
%                 plt = sig(idx:idxe);
%                 pltH = max(plt); pltL = min(plt);
%                 pltRH = max(pltR); pltRL = min(pltR);
%                 scl = (pltH - pltL)/(pltRH - pltRL);
%                 ofst = pltL - pltRL;
            end
        end % alignPilot
    end % methods
end