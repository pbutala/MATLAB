%% FUNCTION [PHI PSI D] = getLOSParams(Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,varargin)
% This function returns the Line of Sight PHI and PSI angles  and distance D
% between specified transmitters and receivers
% 
% -INPUT-
% Xtx Ytx Ztx: X Y Z coordinates of the transmitter(s) 
% Xrx Yrx Zrx: X Y Z coordinates of the receiver(s)
%
% Specify input arguments in any order in the specified format
% 1. ElevationTx: elevation of transmitter.
% 2. AzimuthTx: azimuth of transmitter.
% 3. ElevationRx: elevation of receiver.
% 4. AzimuthRx: azimuth of receiver.
% 
% -OUTUT-
% PHI: Angle subtended between the tx-rx LOS vector and normal to
% transmitter
% PSI: Angle subtended between the rx-tx LOS vector and normal to
% receiver
% D: Length of the tx-rx LOS vector 
% 
% TODO: improve speed if Azi and Ele are 0s (maybe two different sub
% functions ?)

function [PHI,PSI,D] = getLOSParams(Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,varargin)
% default initialization for varargin params

ZenTx = 0;      % zenith Tx
AziTx = 0;      % azimuth Tx
TiltTx = 0;     % tilt Tx
ZenRx = 0;      % zenith Rx
AziRx = 0;      % azimuth Rx
TiltRx = 0;     % tilt Rx

% varargin check
nVArg = nargin-6;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    
    switch lower(varargin{ArgName})
        case 'zenithtx'
            if(isscalar(varargin{ArgParam}))
                ZenTx = varargin{ArgParam};
            else
                error('ZenithTx must be scalar');
            end
        case 'azimuthtx'
            if(isscalar(varargin{ArgParam}))
                AziTx = varargin{ArgParam};
            else
                error('AzimuthTx must be scalar');
            end
        case 'tilttx'
            if(isscalar(varargin{ArgParam}))
                TiltTx = varargin{ArgParam};
            else
                error('TiltTx must be scalar');
            end
        case 'zenithrx'
            if(isscalar(varargin{ArgParam}))
                ZenRx = varargin{ArgParam};
            else
                error('ZenithRx must be scalar');
            end
        case 'azimuthrx'
            if(isscalar(varargin{ArgParam}))
                AziRx = varargin{ArgParam};
            else
                error('AzimuthRx must be scalar');
            end
        case 'tiltrx'
            if(isscalar(varargin{ArgParam}))
                TiltRx = varargin{ArgParam};
            else
                error('TiltRx must be scalar');
            end
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(~isvector(Xtx))
    error('Xtx Ytx Ztx Xrx Yrx Zrx must be vectors');
end

if(~isequal(size(Xtx),size(Ytx),size(Ztx)))
    error('Xtx Ytx Ztx dimensions must agree');
end

if(~isequal(size(Xrx),size(Yrx),size(Zrx)))
    error('Xrx Yrx Zrx dimensions must agree');
end

%% Logic
Ntx = numel(Xtx);
Nrx = numel(Xrx);

PHI = zeros(Nrx,Ntx);
PSI = zeros(Nrx,Ntx);
D = zeros(Nrx,Ntx);
coTr2g = eye(3);

% Nt = [0;0;-1];
Nt = [0;0;1];
% Ntea = getVectRot(Nt,eleTx,'elevation');
% Ntea = getVectRot(Ntea,aziTx,'azimuth');
Ntea = getRotation('Global',Nt,ZenTx,AziTx,TiltTx,coTr2g);

Nr = [0;0;1];
% Nrea = getVectRot(Nr,eleRx,'elevation');
% Nrea = getVectRot(Nrea,aziRx,'azimuth');
Nrea = getRotation('Global',Nr,ZenRx,AziRx,TiltRx,coTr2g);
irx = 1;
while(irx <= Nrx)
    Xrxi = Xrx(irx);
    Yrxi = Yrx(irx);
    Zrxi = Zrx(irx);
    
    itx = 1;
    while(itx <= Ntx)    
        % Initialize each case with specified parameters
        sXtr = Xrxi-Xtx(itx);
        sYtr = Yrxi-Ytx(itx);
        sZtr = Zrxi-Ztx(itx);
        sVtr = [sXtr;sYtr;sZtr];
        
        PHI(irx,itx) = getVtxAng(sVtr,Ntea);
        PSI(irx,itx) = getVtxAng(-sVtr,Nrea);
        D(irx,itx) = norm(sVtr);
        
        itx = itx + 1;
    end
    irx = irx + 1;
end
clearvars -except PHI PSI D;

