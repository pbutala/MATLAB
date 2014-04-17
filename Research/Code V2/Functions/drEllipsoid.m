% drEllipsoid(H,Xloc,Yloc,Zloc,Xrad,Yrad,Zrad,varargin)
% Draws an ellipsoid using the 'surf' function. 
% 
% H: Handle of figure to draw the ellipsoid on.
%
% Xloc Yloc Zloc: X Y Z co-ordinates of the center of the ellipsoids.
%
% Xrad Yrad Zrad: X Y Z radii of the ellipsoids.
%
% 'Color': Sets the color of the ellipsoid. Must be a 1x3 RGB array
%
% 'Alpha': Sets the alpha of the ellipsoid (0-1) with 1=opaque
%
function drEllipsoid(Xloc,Yloc,Zloc,Xrad,Yrad,Zrad,varargin)
% default initialization for varargin params
fClr = [0 0 1];
fAlpha = 1;
zeta = 0;
alpha = 0;
tau = 0;
coTr2g = eye(3);
% Xro = 0;
% Yro = 0;
% Zro = 0;

% varargin check
nVArg = nargin-6;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    switch lower(varargin{ArgName})
        case 'zenith'
            zeta = varargin{ArgParam};
        case 'azimuth'
            alpha = varargin{ArgParam};
        case 'tilt'
            tau = varargin{ArgParam};
        case 'color'
            fClr = varargin{ArgParam};
        case 'alpha'
            fAlpha = varargin{ArgParam};
%         case 'xrotoff'
%             Xro = varargin{ArgParam};
%         case 'yrotoff'
%             Yro = varargin{ArgParam};
%         case 'zrotoff'
%             Zro = varargin{ArgParam};
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(numel(Xrad)==1)
    tXrad = repmat(Xrad,size(Xloc));
    tYrad = repmat(Yrad,size(Yloc));
    tZrad = repmat(Zrad,size(Zloc));
else
    tXrad = Xrad;
    tYrad = Yrad;
    tZrad = Zrad;
end

if(numel(zeta)==1)
    zeta = repmat(zeta,size(Xloc));
end

if(numel(alpha)==1)
    alpha = repmat(alpha,size(Xloc));
end

if(numel(tau)==1)
    tau = repmat(tau,size(Xloc));
end

% if(numel(Xro)==1)
%     Xro = repmat(Xro,size(Xloc));
% end
% 
% if(numel(Yro)==1)
%     Yro = repmat(Yro,size(Xloc));
% end
% 
% if(numel(Zro)==1)
%     Zro = repmat(Zro,size(Xloc));
% end

Ne = numel(Xloc);
idx = 1;

% set(0,'CurrentFigure',H);
hold all;
while(idx <= Ne)
    [X,Y,Z] = ellipsoid(Xloc(idx),Yloc(idx),Zloc(idx),...
                        tXrad(idx),tYrad(idx),tZrad(idx));
    % compute elevation
%     [Xr Yr Zr] = getVectRot(X-Xloc,Y-Yloc,Z-Zloc,ele(idx),'elevation');
    [Xr Yr Zr] = getRotation('Global',X(:)',Y(:)',Z(:)',zeta(idx),alpha(idx),tau(idx),coTr2g,Xloc(idx),Yloc(idx),Zloc(idx),false);
    % compute azimuth
%     [Xr Yr Zr] = getVectRot(Xr,Yr,Zr,azi(idx),'azimuth');
    
    % get X Y Z
%     Xr = Xr+Xloc;
%     Yr = Yr+Yloc;
%     Zr = Zr+Zloc;
    % draw
    Xr = reshape(Xr,size(X));
    Yr = reshape(Yr,size(Y));
    Zr = reshape(Zr,size(Z));
    surf(Xr,Yr,Zr,'EdgeColor','none','FaceColor',fClr,'FaceAlpha',fAlpha);
    
    idx = idx + 1;
end
hold off;
