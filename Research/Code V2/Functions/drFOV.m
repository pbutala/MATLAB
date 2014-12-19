%% Comments
% drFOV(H,Xrx,Yrx,Zrx,Vh,Xlmin,Xlmax,Ylmin,Ylmax,Zlmin,Zlmax,varargin)
% Draws a Cone using the 'patch' function.
%
% H: Handle of figure to draw the Cone on.
%
% Xrx Yrx Zrx: X Y Z co-ordinates of the center of the receivers
% Xl Yl Zl specify limits to draw cone within.
%
% Specify optional input arguments in any order.
%
% 1.'ColorSide': Akin to second/second-fifth row(s) of 'color'
%
% 2.'Elevation': Rotation of the Cone about the X-axis. default = 0rad
% wrt +Y-Axis
%
% 3.'Azimuth': Rotation of the Cone about the Z-axis. default = 0rad wrt
% +X-axis Rotation logic implements elevation first, then azimuth. Rotation
% about Y-axis can be implemented by a combination of the other two.
%
% 4.'SideAlpha': Sets the alpha of the sides (0-1) with 1=opaque
%
%% Function
function drFOV(Xrx,Yrx,Zrx,Vh,Xlmin,Xlmax,Ylmin,Ylmax,Zlmin,Zlmax,varargin)
% default initialization for varargin params
Nsides = 90;
FclrS = [1 1 0];
sAlpha = 0.5;
zeta = 0;
alpha = 0;
tau = 0;

% varargin check
nVArg = nargin-10;
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
        case 'colorside'
            FclrS = varargin{ArgParam};
        case 'sidealpha'
            sAlpha = varargin{ArgParam};
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if ~(isvector(FclrS) && (numel(FclrS) == 3))
    error('ColorSide should be a 3 element vector');
end

if(numel(zeta)==1)
    zeta = repmat(zeta,size(Xrx));
end

if(numel(alpha)==1)
    alpha = repmat(alpha,size(Xrx));
end

if(numel(sAlpha)==1)
    sAlpha = repmat(sAlpha,size(Xrx));
end

if(isscalar(Vh))
    Vh = Vh*ones(size(Xrx));
end
if(isscalar(Xlmin))
    Xlmin = Xlmin*ones(size(Xrx));
    Xlmax = Xlmax*ones(size(Xrx));
    Ylmin = Ylmin*ones(size(Xrx));
    Ylmax = Ylmax*ones(size(Xrx));
    Zlmin = Zlmin*ones(size(Xrx));
    Zlmax = Zlmax*ones(size(Xrx));
end
if(~isequal(size(Xlmin),size(Ylmin),size(Zlmin),...
        size(Xlmax),size(Ylmax),size(Zlmax),size(Xrx)))
    error('Xl Yl Zl min max dimensions must agree with Xrx');
end

%% logic
NXc = 2*(Nsides+1);
FcsS = [2:2:NXc-2; 3:2:NXc-1; 4:2:NXc]';

Rlim = sqrt((Xlmin-Xlmax).^2 +...
    (Ylmin-Ylmax).^2 +...
    (Zlmin-Zlmax).^2);

Np = numel(Xrx);
idx = 1;

% set(0,'CurrentFigure',H);

while(idx <= Np)
    Xmin = Xlmin(idx);
    Xmax = Xlmax(idx);
    Ymin = Ylmin(idx);
    Ymax = Ylmax(idx);
    Zmin = Zlmin(idx);
    Zmax = Zlmax(idx);
    Rmax = Rlim(idx);
    
    
    [Xc,Yc,Zc] = cylinder([0 tan(Vh(idx))],Nsides);
    Xc = Xc*Rmax;
    Yc = Yc*Rmax;
    Zc = Zc*Rmax;
    
    Vtx = [Xc(:) Yc(:) Zc(:)];
    
%     % compute elevation
%     Vtx = getVectRot(Vtx',zeta(idx),'elevation')';
%     
% %     % compute azimuth
%     Vtx = getVectRot(Vtx',alpha(idx),'azimuth')';
    Vtx = getRotation('Global',Vtx',zeta(idx),alpha(idx),tau(idx),eye(3));
    % compute offsets
    ofst = repmat([Xrx(idx) Yrx(idx) Zrx(idx)],[size(Vtx',1) 1]);
    Vtx = Vtx' + ofst;
    
    % compute min max constraints.
    % even rows are base vertex
    for r=2:2:size(Vtx,1)
        v = Vtx(r,:) - Vtx(r-1,:);
        u = v/norm(v);
        
        scX1 = (Xmin-Xrx(idx))/u(1);
        scX2 = (Xmax-Xrx(idx))/u(1);
        scY1 = (Ymin-Yrx(idx))/u(2);
        scY2 = (Ymax-Yrx(idx))/u(2);
        scZ1 = (Zmin-Zrx(idx))/u(3);
        scZ2 = (Zmax-Zrx(idx))/u(3);
        
        srtSc = sort([scX1 scX2 scY1 scY2 scZ1 scZ2],2);
        sc = srtSc(find(srtSc>0, 1));
        
        Vtx(r,:) = Vtx(r-1,:) + sc*u;
    end
    
    patch('Vertices',Vtx, 'Faces', FcsS,...
        'EdgeColor', FclrS,'EdgeAlpha',sAlpha(idx),...
        'FaceColor', FclrS,'FaceAlpha',sAlpha(idx));
    
    idx = idx + 1;
end

