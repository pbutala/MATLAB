classdef cDemodOFDM < cDemodulator
    % class to handle Optical OFDM demodulation
    properties(Constant)
        CLIPLOW = 0;
        CLIPHIGH = realmax('double');
    end % properties - constant
    
    properties(SetAccess = immutable)
        OFDMTYP;                    % Type of optical OFDM
        NSC;                        % Number of subcarriers
        SYMSC;                      % Subcarrier Symbols (Typically QAM)
        MSC;                        % Modulation order for subcarrier modulation
        OFST;                       % DC offset
    end % properties - immutable
    
    properties(SetAccess = protected)
        BPSYM;                      % Bits Per Symbol
        NPSYM;                      % Samples Per Symbol
    end % properties - protected
    
    properties
        FILTER = 'RAISEDCOSINE';       % Up/Dn sampling filter type
    end % properties
    
    methods
        % CONSTRUCTOR
        function obj = cDemodOFDM(ofdmtyp, nsc, syms, ofst,...
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
            nps = ceil(nsc*(clkin/clkout)); % samples per symbol
            if~exist('clkin','var')
                clkin = 1;                  % default input clock
            end
            if~exist('clkout','var')
                clkout = 1;                 % default output clock
            end
            if~exist('bufszin','var')
                bufszin = nps;             % default input buffer size
            end
            if~exist('bufszout','var')
                bufszout = bps;      % default output buffer size
            end
            obj = obj@cDemodulator(clkin, clkout, bufszin, bufszout);
            obj.OFDMTYP = ofdmtyp;
            obj.NSC = nsc;
            obj.SYMSC = syms;
            obj.MSC = msc;
            obj.OFST = ofst;
            obj.BPSYM = bps;
            obj.NPSYM = nps;
        end % constructor
    end % methods
    
    methods 
        % MODULATE
        function demodulate(obj)
            while(obj.BUFIN.COUNT >= obj.NPSYM)
                strm = obj.BUFIN.deQ(obj.NPSYM);
                tSig = updnClock(strm,obj.CLKIN,obj.CLKOUT,obj.FILTER,false);
                rxDat = decodeOFDMsignal(tSig,...
                    'OFDMtype',obj.OFDMTYP,...
                    'N',obj.NSC,...
                    'Symbols',obj.SYMSC,...
                    'ShowRcv',true);
                tBits = dec2binMat(rxDat-1, log2(obj.MSC));
                tBits = tBits.';
                obj.BUFOUT.enQ(tBits(:));
            end % while loop
        end % demodulate
    end % methods - overloaded
end






