classdef cModOFDM < cModulator
    % class to handle Optical OFDM modulation
    
    properties(SetAccess = immutable)
        OFDMTYP;                    % Type of optical OFDM
        NSC;                        % Number of subcarriers
        SYMSC;                      % Subcarrier Symbols (Typically QAM)
        MSC;                        % Modulation order for subcarrier modulation
        DCOOFST;                    % DC offset
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
        function obj = cModOFDM(ofdmtyp, nsc, syms, dcoofst, clpLo, clpHi, scl, ofst,...
                clkin, clkout, bufszin, bufszout)
            msc = numel(syms);              % subcarrier modulation order (M-QAM)
            switch lower(ofdmtyp)
                case 'acoofdm'
                    d = nsc/4;                % number of data carriers per ACOOFDM symbol
                case {'dcoofdm','dmt'}
                    d = nsc/2 - 1;            % number of data carriers per DCOOFDM symbol
                otherwise
                    error('OFDM type must be ''ACOOFDM'' or ''DCOOFDM'' or ''DMT''');
            end
            bps = log2(msc)*d;              % bits per symbol
            nps = ceil(nsc*clkout/clkin);   % samples per symbol
            if~exist('clpLo','var')
                clpLo = 0;                  % default clip lo level
            end
            if~exist('clpHi','var')
                clpHi = realmax('double');  % default clip hi level
            end
            if~exist('scl','var')
                scl = 1;
            end
            if~exist('ofst','var')
                ofst = 0;
            end
            if~exist('clkin','var')
                clkin = 1;                  % default input clock
            end
            if~exist('clkout','var')
                clkout = 1;                 % default output clock
            end
            if~exist('bufszin','var')
                bufszin = bps;              % default input buffer size
            end
            if~exist('bufszout','var')
                bufszout = nps;             % default output buffer size
            end
            obj = obj@cModulator(clkin, clkout, bufszin, bufszout, clpLo, clpHi, scl, ofst);
            obj.OFDMTYP = ofdmtyp;
            obj.NSC = nsc;
            obj.SYMSC = syms;
            obj.MSC = msc;
            obj.DCOOFST = dcoofst;
            obj.BPSYM = bps;
            obj.NPSYM = nps;
        end % constructor
    end % methods
    
    methods 
        % MODULATE
        function sig = modulate(obj)
            while(obj.BUFIN.COUNT >= obj.BPSYM)
                strm = obj.BUFIN.deQ(obj.BPSYM);
                tBits = reshape(strm,[log2(obj.MSC), obj.BPSYM/log2(obj.MSC)]);
                data = bin2decMat(tBits.') + 1;
                % Generate OFDM symbol
                tSig = genOFDMsignal(... % Variable Arguments to the function
                    'data',data,...
                    'OFDMtype',obj.OFDMTYP,...
                    'N',obj.NSC,...
                    'Symbols',obj.SYMSC,...
                    'ClipLow',obj.SIGLO,...
                    'ClipHigh',obj.SIGHI,...
                    'OffsetDcoStddev', obj.DCOOFST,...
                    'ShowConst',false);
                sig = updnClock(tSig,obj.CLKIN,obj.CLKOUT,obj.FILTER,false);
            end % while loop
        end % modulate
    end % methods - overloaded
end






