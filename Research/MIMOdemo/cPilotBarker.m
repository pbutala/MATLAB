classdef cPilotBarker < cPilot
    % base class to handle Barker pilots.
    properties(SetAccess = immutable)
        TYPE;
    end % properties - immutable
    
    properties
        FILTER = 'RAISEDCOSINE';    % Up/Dn sampling filter type
    end % properties
    
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
                case 'BARKER1'
                    plt = 1;
                case 'BARKER2'
                    plt = [1;0];
                case 'BARKER3'
                    plt = [1;1;0];
                case 'BARKER4'
                    plt = [1;1;0;1];
                case 'BARKER5'
                    plt = [1;1;1;0;1];
                case 'BARKER7'
                    plt = [1;1;1;0;0;1;0];
                case 'BARKER11'
                    plt = [1;1;1;0;0;0;1;0;0;1;0];
                case 'BARKER13'
                    plt = [1;1;1;1;1;0;0;1;1;0;1;0;1];
            end
            val = updnClock(plt, obj.CLKIN, clkout, obj.FILTER, false);
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
%                 figure;
%                 plot(lag,acor);
%                 drCrosshair(lag(I), acor(I), true, 'half', 'r');
%                 str = sprintf(' lag = %d',lag(I));
%                 text(lag(I), acor(I), str);
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
            idx = IDXs(find(MSE==max(MSE),1,'first'));
        end % alignFine
    end % methods - private
end



























































