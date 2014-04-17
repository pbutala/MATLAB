%% Comments
% drPyramid(H,Xloc,Yloc,Zloc,varargin) 
% Draws a pyramid using the 'patch' function.
% Default: Pyramid has base in the X-Y plane and apex along
% +Z-axis. Each side of base and height is 1 unit.
%
% H: Handle of figure to draw the pyramid on.
%
% Xloc Yloc Zloc: X Y Z co-ordinates of the center of the bases to position
% the pyramids.
%
% Specify optional input arguments in any order.
%
% 1.'Scale': scaling factor for the pyramid. Can provide one scale for all
% or a vector of scales, one for each pyramid.
%
% 2.'Direction': + values indicate upward pointing pyramid while -ve values
% indicate downward pointing pyramid. The magnitude of the value scales the
% height wrt the edge.
%
% 3.'Color': Can be a 2x3 or 5x3 vector per pyramid.For multiple pyramids
% of different colors, an array of vectors can be specified. First row
% specifies the RGB values for color of the base. For 2x3, second row
% specifies common color for all faces. For 5x3, second-fifth rows specify
% color for faces spanning 1-2,2-3,3-4,4-1 quadrants at default. Default
% colors are black for base and yellow for sides.
%
% 4.'ColorBase': Akin to first row of 'color'
%
% 5.'ColorSide': Akin to second/second-fifth row(s) of 'color'
%
% 6.'Elevation': Rotation of the pyramid about the X-axis. default = 0rad
% wrt +Y-Axis
%
% 7.'Azimuth': Rotation of the pyramid about the Z-axis. default = 0rad wrt
% +X-axis Rotation logic implements elevation first, then azimuth. Rotation
% about Y-axis can be implemented by a combination of the other two.
%
% 8.'BaseAlpha': Sets the alpha of the base (0-1) with 1=opaque
%
% 9.'SideAlpha': Sets the alpha of the sides (0-1) with 1=opaque
%
%% Function
function drPyramid(Xloc,Yloc,Zloc,varargin)
% default initialization for varargin params
scl = 1;
dir = 1;
apex = [];
FclrB = [0 0 0];
FclrS = [1 1 0; 1 1 0; 1 1 0; 1 1 0];
bAlpha = 1;
sAlpha = 1;
zeta = 0;
alpha = 0;
tau = 0;

% varargin check
nVArg = nargin-3;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    switch lower(varargin{ArgName})
        case 'apex'
            apex = varargin{ArgParam};
        case 'scale'
            scl = varargin{ArgParam};
        case 'direction'
            dir = varargin{ArgParam};
        case 'zenith'
            zeta = varargin{ArgParam};
        case 'azimuth'
            alpha = varargin{ArgParam};
        case 'tilt'
            tau = varargin{ArgParam};
        case 'colorbase'
            FclrB = varargin{ArgParam};
        case 'colorside'
            FclrS = varargin{ArgParam};
        case 'color'
            FclrB = varargin{ArgParam}(1,:);
            FclrS = varargin{ArgParam}(2:end,:);
        case 'basealpha'
            bAlpha = varargin{ArgParam};
        case 'sidealpha'
            sAlpha = varargin{ArgParam};
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(numel(scl)==1)
    tScl = repmat(scl,size(Xloc));
else
    tScl = scl;
end

if(numel(dir)==1)
    tDir = repmat(dir,size(Xloc));
else
    tDir = dir;
end

if((numel(FclrB)==3) && (size(FclrB,1)==1))
    tFclrB = repmat(FclrB,[1 numel(Xloc)]);
else
    tFclrB = FclrB;
end

if((numel(FclrS)==12) && (size(FclrS,1)==4))
    tFclrS = repmat(FclrS,[1 numel(Xloc)]);
else
    tFclrS = FclrS;
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

if(numel(bAlpha)==1)
    tbAlpha = repmat(bAlpha,size(Xloc));
else
    tbAlpha = bAlpha;
end

if(numel(sAlpha)==1)
    tsAlpha = repmat(sAlpha,size(Xloc));
else
    tsAlpha = sAlpha;
end

if(~isempty(apex))
if ~(isrow(apex) && (numel(apex)==3))
    error('Apex must be a 1x3 vector');
end
end

%% logic
FcsB = [1 2 3 4];
FcsS = [1 2 3; 1 3 4; 1 4 5; 1 5 2];

Np = numel(Xloc);
idx = 1;
clrR = 1;

% set(0,'CurrentFigure',H);

while(idx <= Np)
    Vtx = [0 0 tDir(idx); 0.5 0.5 0; -0.5 0.5 0; -0.5 -0.5 0; 0.5 -0.5 0];
    
%     % compute elevation
%     Vtx = getVectRot(Vtx',tEle(idx),'elevation')';
% 
%     % compute azimuth
%     Vtx = getVectRot(Vtx',tAzi(idx),'azimuth')';
    Vtx = getRotation('global',Vtx',zeta(idx),alpha(idx),tau(idx),eye(3));
    
    % compute offsets
    ofstB = repmat([Xloc(idx) Yloc(idx) Zloc(idx)],[4,1]);
%     VtxB = tScl(idx)*Vtx(2:end,:)+ofstB;
    VtxB = tScl(idx)*Vtx(:,2:end)'+ofstB;
    
    ofstS = repmat([Xloc(idx) Yloc(idx) Zloc(idx)],[5,1]);
    VtxS = tScl(idx)*Vtx' +ofstS;
    % check Apex
    if(~isempty(apex))
        VtxS(1,:) = apex;
    end
    
    % coloring is weird but thats how this works
    patch('Vertices',VtxB, 'Faces', FcsB,...
         'FaceColor', tFclrB(clrR:clrR+2)', 'FaceAlpha',tbAlpha(idx));
    patch('Vertices',VtxS, 'Faces', FcsS,...
         'FaceVertexCData', tFclrS(:,clrR:clrR+2),'FaceColor', 'flat',...
         'FaceAlpha',tsAlpha(idx));
    
    clrR = clrR + 3;
    idx = idx + 1;
end
