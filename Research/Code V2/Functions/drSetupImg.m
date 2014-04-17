%% Comment
% drSetupImg(H,aXrm,aYrm,aZrm,Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,varargin)
% Specify optional input arguments in any order.
% 1. ColorBaseTx: color of base of transmitters (default [0 0 0] black)
% 2. FOVhalf: Field of View half angle (default pi/4)
% 3. LightRays: 'None' 'FOV' and 'Pyramid' are options that control how the
% light rays from the transmitter to receiver are shown. (dafault FOV)

%% Function
function drSetupImg(aXrm,aYrm,aZrm,Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,varargin)
% default initialization for varargin params
clrBtx = [0 0 0];
zetaTx = 0;
alphaTx = 0;
tauTx = 0;
zetaRx = 0;
alphaRx = 0;
tauRx = 0;
Vh = pi/4;
lRays = 0;
% varargin check
nVArg = nargin-9;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    switch lower(varargin{ArgName})
        case 'colorbasetx'
            clrBtx = varargin{ArgParam};
        case 'zenithtx'
            zetaTx = varargin{ArgParam};
        case 'azimuthtx'
            alphaTx = varargin{ArgParam};
        case 'tilttx'
            tauTx = varargin{ArgParam};
        case 'zenithrx'
            zetaRx = varargin{ArgParam};
        case 'azimuthrx'
            alphaRx = varargin{ArgParam};
        case 'tiltrx'
            tauRx = varargin{ArgParam};
        case 'fovhalf'
            Vh = varargin{ArgParam};
        case 'lightrays'
            switch lower(varargin{ArgParam})
                case 'fov'
                    lRays = 0;
                case 'pyramid'
                    lRays = 1;
                case 'none'
                    lRays = 2;
                otherwise
                     error('unknown parameter %s specified',varargin{ArgParam});
            end
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(~isscalar(aXrm))
    error('aXrm,aYrm,aZrm must be scalar.');
end

if(~isequal(size(aXrm),size(aYrm),size(aZrm)))
    error('aXrm Yrx aZrm dimensions must agree');
end

if(~isequal(size(Xtx),size(Ytx),size(Ztx)))
    error('Xtx Ytx Ztx dimensions must agree');
end

if(~isequal(size(Xrx),size(Yrx),size(Zrx)))
    error('Xrx Yrx Zrx dimensions must agree');
end

if((numel(clrBtx)==3) && (size(clrBtx,1)==1))
    clrBtx = repmat(clrBtx,[1 numel(Xtx)]);
end

% default configuration
atx = 20e-2;
arx = 20e-2;% just so it is seen on the figure
Do = 20e-2; % just so it is seen on the figure
to = 1e-2;% just so it is seen on the figure
f = 1e-1; % just so it is seen on the figure

%% Logic
% set(0,'CurrentFigure',H);
axis([0 aXrm 0 aYrm 0 aZrm]);
axis equal;
% Draw transmitters
drPyramid(Xtx,Ytx,Ztx,'Scale',atx,'direction',1,...
          'colorbase',clrBtx,'zenith',zetaTx,'azimuth',alphaTx,'tilt',tauTx);

% Draw lens
drEllipsoid(Xrx,Yrx,Zrx,Do,Do,to,'alpha',0.25,...
    'zenith',zetaRx,'azimuth',alphaRx,'tilt',tauRx);

% Draw sensor
drPyramid(Xrx,Yrx,Zrx-f,'Scale',arx,'direction',0.25,...
    'basealpha',0.5,'sidealpha',0.5,'zenith',zetaRx,'azimuth',alphaRx,'tilt',tauRx);

% Draw light rays
switch lRays
    case 0
        drFOV(Xrx,Yrx,Zrx,Vh,...
            0, aXrm,0, aYrm,0, max(Ztx(:)),...
        'zenith',zetaRx,'azimuth',alphaRx,'tilt',tauRx);
        
    case 1
        Scale = (max(Xtx(:))-min(Xtx(:)));
        drPyramid(aXrm/2,aYrm/2,min(Ztx(:)-atx),...
            'Apex',[Xrx Yrx Zrx],...
            'Scale',Scale,...
            'colorbase',[1,1,1],'basealpha',0.5,'sidealpha',0.5,...
            'zenith',zetaRx,'azimuth',alphaRx,'tilt',tauRx);
end

% grid on;
% view(3);
% xlabel('X: length');
% ylabel('Y: width');
% zlabel('Z: height');

%% Space






































