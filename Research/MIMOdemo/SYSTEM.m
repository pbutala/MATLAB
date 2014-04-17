classdef SYSTEM
    properties(SetAccess = private)
        % FRAME
        frmPILOTCYCLES = 4;            % pilot cycles
        frmPILOTf = 10e6;             % pilot sin wave frequency (Hz)
        frmPILOTTYPE = 'BARKER11';      % pilot signal type
        frmPILOTSF = 1;                 % for barker code only. pilot sample scale factor (without this, the pilot is assumed to have been sampled at signal sample rate.)
        % OFDM
        dFs = 40e6;                     % Signal: sample rate (2CE the signal BW)
        modTYPE = 'DCOOFDM';            % type of OFDM
        modOFST = 3.2;                  % offset applied to time domain signal
        modN = 128;              % # total subcarriers
        modM = 4;                       % M-QAM
        modCLIPH = 3.2;                   % clip high
        modCLIPL = 0;                   % clip low
        
        % OTHER
        ctPI = pi;
        cti = sqrt(-1);
    end
    
    properties(Dependent = true, SetAccess = private)
        % FRAME
        frmPILOTLEN16;            % pilot length in int16
        frmPILOTLEN16SF;            % length of up-/down- scaled PILOT signal in number of ints
        frmLEN16;                       % length of frame in number of unsigned ints
        frmLEN16SF;                 % length of up-/down- scaled frame in number of unsigned ints
        frmLEN8;                        % length of frame in number of bytes (unsigned chars)
        frmLEN8SF;                        % length of up-/down- scaled frame in number of bytes (unsigned chars)
        frmDATALEN16;                   % # data int16 in each frame
        frmPILOT01;                     % PILOT signal (min:0, max:1) based on properties specified
        % OFDM
        % DEVICE
        dSF;                            % device up/down - scale factor.
    end
    
    properties(Abstract, Dependent = true, SetAccess = private)
        % DEVICE
        dETHID; % device network card id
        dCLKs;                   % device transmit/receive sample clock
        dCLKSRC;                   % device clock source for interface
    end
    
    methods(Static, Access=protected)
        function val = GetNetPCIcards(strType)
            [~,str] = system('"C:\Program Files (x86)\Windows Kits\8.1\Tools\x64\devcon" find =Net PCI\*');
            iNA = strfind(str,'4&11050A08&0&00E5');
            iTX = strfind(str,'4&18BDCD5B&0&48F0');
            iRX = strfind(str,'4&18BDCD5B&0&20F0');
            idxs = sort([iNA;iTX;iRX],1,'descend');
            switch(lower(strType))
                case 'tx'
                val = find(idxs == iTX)-1;
                case 'rx'
                val = find(idxs == iRX)-1;
            end
        end
    end
    % FUNCTIONS
    methods
        % returns the PILOT signal (min:0, max:1) based on properties specified
        % DEFAULT: SINE
        function val = getPILOT01(obj)
            plt = sin((0:obj.frmPILOTLEN16-1)'*2*obj.ctPI*obj.frmPILOTf/obj.dFs);
            switch(lower(obj.frmPILOTTYPE))
                case 'sine'
                    plt = (plt+1)*0.5;
                case 'square'
                    plt(plt>=0.5) = 1;
                    plt(plt<0.5) = 0;
                case 'barker11'
                    plt = [1;1;1;0;0;0;1;0;0;1;0];
                    pltR = repmat(plt',obj.frmPILOTSF,1);
                    plt = pltR(:);
                case 'barker13'
                    plt = [1;1;1;1;1;0;0;1;1;0;1;0;1];
                    pltR = repmat(plt',obj.frmPILOTSF,1);
                    plt = pltR(:);
                otherwise
                    error('Pilot type must be ''SINE'',''SQUARE'',''BARKER11'',''BARKER13''');
            end
            val = plt;
            % filter pilot to avoid aliasing in data signal
%             [Fb,Fa] = butter(4,0.8);
%             val = filter(Fb,Fa,plt);
        end
    end
    % GETTERS/SETTERS
    methods
        % returns length of PILOT signal in number of ints
        function val = get.frmPILOTLEN16(obj)
            switch(lower(obj.frmPILOTTYPE))
                case {'sine','square'} 
                    val = ceil(obj.dFs*obj.frmPILOTCYCLES/obj.frmPILOTf);
                case {'barker11'}
                    val = 11*obj.frmPILOTSF;
                case {'barker13'}
                    val = 13*obj.frmPILOTSF;
                otherwise
                    error('Pilot type must be ''SINE'',''SQUARE'',''BARKER11'',''BARKER13''');
            end
        end
        
        % returns length of up-/down- scaled PILOT signal in number of ints
        function val = get.frmPILOTLEN16SF(obj)
            val = ceil(obj.frmPILOTLEN16*obj.dSF);
        end
        % returns TOTAL length of frame in number of unsigned ints
        function val = get.frmLEN16(obj)
            val = obj.frmPILOTLEN16 + obj.modN;
        end
        
        % returns TOTAL length of up-/down- scaled frame in number of unsigned ints
        function val = get.frmLEN16SF(obj)
            val = ceil(obj.frmLEN16*obj.dSF);
        end
        
        % returns TOTAL length of frame in number of bytes (unsigned chars)
        function val = get.frmLEN8(obj)
            val = obj.frmLEN16*2;
        end
        
        % returns TOTAL length of up-/down- scaled frame in number of bytes (unsigned chars)
        function val = get.frmLEN8SF(obj)
            val = obj.frmLEN16SF*2;
        end
        
        % gets number of data samples in each frame
        function val = get.frmDATALEN16(obj)
            switch lower(obj.modTYPE)
                case 'acoofdm'
                    val = obj.modN/4;    % number of data carriers per ACOOFDM symbol
                case {'dcoofdm','dmt'}
                    val = obj.modN/2 - 1;    % number of data carriers per DCOOFDM symbol
                otherwise
                    error('OFDM type must be ''ACOOFDM'' or ''DCOOFDM'' or ''DMT''');
            end
        end
        
        % gets transmitter UpSample Factor
        function val = get.dSF(obj)
            val = obj.dCLKs/obj.dFs;
        end
    end
end







































