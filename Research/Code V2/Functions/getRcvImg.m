%% COMMENTS
% [Xs1 Ys1 Xs2 Ys2 Xs3 Ys3 Xs4 Ys4] = getRcvImg(Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,atx,f,varargin)
% Calculates the center location and the side length of spots on image.
% Also returns if the transmitter corresponding to a spot is in FOV and at
% what angle psi.
%
% If the transmitter corresponding to a spot is outside the receiver's FOV,
% the function will return the location [Xs Ys Zs] as if the receivers FOV
% half was pi/2 BUT the spot side length as will be set to 0.
%
% Specify input arguments in any order in the specified format
% 1. FOVHalf: Vield of View half angle. (default pi/2)

%% Function
function [Xs1,Ys1,Xs2,Ys2,Xs3,Ys3,Xs4,Ys4] = getRcvImg(Xtx,Ytx,Ztx,Xrx,Yrx,Zrx,atx,f,varargin)
% default initialization for varargin params
Vh = pi/2;
zetaTx = 0;
alphaTx = 0;
tauTx = 0;
zetaRx = 0;
alphaRx = 0;
tauRx = 0;

% varargin check
nVArg = nargin-8;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    switch lower(varargin{ArgName})
        case 'fovhalf'
            Vh = varargin{ArgParam};
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
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(~isequal(size(Xtx),size(Ytx),size(Ztx)))
    error('Xtx Ytx Ztx dimensions must agree');
end

if(~isequal(size(Xrx),size(Yrx),size(Zrx)))
    error('Xrx Yrx Zrx dimensions must agree');
end

if (isscalar(atx))
    atx = atx*ones(size(Xtx));
elseif (isequal(size(atx),size(Xtx)))
    %     atx = atx;
else
    error('Transmitter side length should be a scalar or matrix of size(Xrx)');
end

if (isscalar(f))
    f = f*ones(size(Xrx));
elseif (isequal(size(f),size(Xrx)))
    %     f = f;
else
    error('Focus should be a scalar or matrix of size(Xrx)');
end

if(numel(Vh)==1)
    Vh = repmat(Vh,size(Xrx));
elseif(~isequal(size(Vh),size(Xrx)))
    error('FOVhalf must be a scalar or a matrix with same dimensions as Xrx');
end

if(numel(zetaTx)==1)
    zetaTx = repmat(zetaTx,size(Xtx));
elseif(~isequal(size(zetaTx),size(Xtx)))
    error('ZenithTx must be a scalar or a matrix with same dimensions as Xtx');
end

if(numel(alphaTx)==1)
    alphaTx = repmat(alphaTx,size(Xtx));
elseif(~isequal(size(alphaTx),size(Xtx)))
    error('AzimuthTx must be a scalar or a matrix with same dimensions as Xtx');
end

if(numel(tauTx)==1)
    tauTx = repmat(tauTx,size(Xtx));
elseif(~isequal(size(tauTx),size(Xtx)))
    error('TiltTx must be a scalar or a matrix with same dimensions as Xtx');
end

if(numel(zetaRx)==1)
    zetaRx = repmat(zetaRx,size(Xrx));
elseif(~isequal(size(zetaRx),size(Xrx)))
    error('ZenithRx must be a scalar or a matrix with same dimensions as Xrx');
end

if(numel(alphaRx)==1)
    alphaRx = repmat(alphaRx,size(Xrx));
elseif(~isequal(size(alphaRx),size(Xrx)))
    error('AzimuthRx must be a scalar or a matrix with same dimensions as Xrx');
end

if(numel(tauRx)==1)
    tauRx = repmat(tauRx,size(Xrx));
elseif(~isequal(size(tauRx),size(Xrx)))
    error('TiltRx must be a scalar or a matrix with same dimensions as Xrx');
end

%% Logic
Ntx = numel(Xtx);
Nrx = numel(Xrx);

% szRxL = size(Xrx,1);
% szRxW = size(Xrx,2);
% szRxH = size(Xrx,3);
% szRx = [szRxL szRxW szRxH];

szTxL = size(Xtx,1);
szTxW = size(Xtx,2);
szTxH = size(Xtx,3);
szTx = [szTxL szTxW szTxH];

% initialize return variables
Xs1 = zeros([Nrx Ntx]);
Ys1 = zeros([Nrx Ntx]);
Xs2 = zeros([Nrx Ntx]);
Ys2 = zeros([Nrx Ntx]);
Xs3 = zeros([Nrx Ntx]);
Ys3 = zeros([Nrx Ntx]);
Xs4 = zeros([Nrx Ntx]);
Ys4 = zeros([Nrx Ntx]);

Nr = [0;0;1];
irx = 1;
while(irx <= Nrx)
    isFOV = false(szTx);
%     Nrea = getVectRot(Nr,zetaRx(irx),'elevation');
%     Nrea = getVectRot(Nrea,alphaRx(irx),'azimuth');
    Nrea = getRotation('global',Nr,zetaRx(irx),alphaRx(irx),tauRx(irx),eye(3));
    itx = 1;
    while(itx <= Ntx)
        % Calculate vector from receiver to transmitter
        Xrt = Xtx(itx) - Xrx(irx);
        Yrt = Ytx(itx) - Yrx(irx);
        Zrt = Ztx(itx) - Zrx(irx);
        Vrt = [Xrt;Yrt;Zrt];
            
            % get co-ordinates of four corners of transmitter wrt the
            % center of the receiver
            dtx = atx(itx)/2;
            Xp = [dtx -dtx -dtx dtx];
            Yp = [dtx dtx -dtx -dtx];
            Zp = [0 0 0 0];
%             [Xp Yp Zp] = getVectRot(Xp,Yp,Zp,zetaTx(itx),'elevation');
%             [Xp Yp Zp] = getVectRot(Xp,Yp,Zp,alphaTx(itx),'azimuth');
            [Xp,Yp,Zp] = getRotation('global',Xp,Yp,Zp,zetaTx(itx),alphaTx(itx),tauTx(itx),eye(3));
            Xp = Xp + Xrt;
            Yp = Yp + Yrt;
            Zp = Zp + Zrt;
            
            % Get co-ordinates of transmitter wrt receiver.
            % this considers the rceiver elevation and azimuth
%             [Xp,Yp,Zp] = getCoord3D(Xp,Yp,Zp,-zetaRx(irx),-alphaRx(irx),-tauRx(irx));
            [Xp,Yp,Zp] = getCoord3D(Xp,Yp,Zp,zetaRx(irx),alphaRx(irx),tauRx(irx));
            
            % Ensure tx and rx are facing each other
            if isequal(sign(Zp),ones(size(Zp)))
                % calculate angle psi.
                psirx = getVtxAng(Vrt,Nrea);
                % if in FOV
                if psirx < Vh
                    isFOV(itx) = true;
                    % Get the coordinates of the image of transmitter in the new
                    % co-ordinate system i.e. origin is at center of receiver
                    [Xs1(irx,itx),Ys1(irx,itx)] = ...
                        getSimplePoint(Xp(1),Yp(1),Zp(1),f(irx));
                    [Xs2(irx,itx),Ys2(irx,itx)] = ...
                        getSimplePoint(Xp(2),Yp(2),Zp(2),f(irx));
                    [Xs3(irx,itx),Ys3(irx,itx)] = ...
                        getSimplePoint(Xp(3),Yp(3),Zp(3),f(irx));
                    [Xs4(irx,itx),Ys4(irx,itx)] = ...
                        getSimplePoint(Xp(4),Yp(4),Zp(4),f(irx));
                end
            end
        itx = itx + 1;
    end
    irx = irx + 1;
end
end

%% Functions

function [Xi,Yi] = getSimplePoint(Xp,Yp,Zp,fcs)
Mag = fcs/(Zp-fcs);
Xi = -Mag*Xp;
Yi = -Mag*Yp;
end
%% Space







































