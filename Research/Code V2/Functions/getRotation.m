% function [Xr Yr Zr coTr2gr] = getRotation('FrOfRef',X,Y,Z,PHI,THETA,PSI,coTr2g,Xp0,Yp0,Zp0,flgRtOfst)
% function [MtXYZr coTr2gr] = getRotation('FrOfRef',MtXYZ,PHI,THETA,PSI,coTr2g,Xp0,Yp0,Zp0,flgRtOfst)
% function [Xr Yr Zr coTr2gr] = getRotation('FrOfRef',X,Y,Z,PHI,THETA,PSI,coTr2g)
% function [MtXYZr coTr2gr] = getRotation('FrOfRef',MtXYZ,PHI,THETA,PSI,coTr2g)
% 
% getRotation function returns XYZ co-ordinates in the global co-ordinate
% system after performing 1.Tilt 2.Elevation 3.Azimuth. This gives same
% results as 1.Elevation 2.Azimuth 3.Tilt 
% 
% -INPUT- 
% FrOfRef: Frame of reference for the rotation co-ordinates. This parameter
% can be either 'Global' or 'Receiver'
% 
% X: X co-ordinates of vectors in global system
% Y: Y co-ordinates of vectors in global system
% Z: Z co-ordinates of vectors in global system
% MtXYZ: Matrix [X;Y;Z]
% 
% PHI: Elevation in radians (scalar)
% THETA: Azimuth in radians (scalar)
% PSI: Tilt in radians (scalar)
%
% coTr2g: Co-ordinate Transformation matrix. Should be 3x3. Each column of
% this matrix should be the receiver basis set as projected on the global
% coordinate system.
% 
% Xp0: X co-ordinates of offset vectors
% Yp0: Y co-ordinates of offset vectors
% Zp0: Z co-ordinates of offset vectors
%
% flgRtOfst: Flag indicating if offset should be rotated.
%            true: rotate offset; false: do not rotate offset
%
% -OUTPUT-
% Xr: X co-ordinates of rotated vectors in global system
% Yr: Y co-ordinates of rotated vectors in global system
% Zr: Z co-ordinates of rotated vectors in global system
% coTr2gr: New co-ordinate transformation matrix.
% 

function varargout = getRotation(varargin)
narginchk(6,12);
Xp0 = 0;
Yp0 = 0;
Zp0 = 0;
flgRtOfst = false;

if(nargin == 6)
% function MtXYZr = getRotation('FrOfRef',MtXYZ,PHI,THETA,PSI,coTr2g)
    frOfRef = varargin{1};
    [X Y Z] = getVtFromMt(varargin{2});
    PHI = varargin{3};
    THETA = varargin{4};
    PSI = varargin{5};
    coTr2g = varargin{6};
elseif(nargin == 8)
% function [Xr Yr Zr] = getRotation('FrOfRef',X,Y,Z,PHI,THETA,PSI,coTr2g)
    frOfRef = varargin{1};    
    X = varargin{2};
    Y = varargin{3};
    Z = varargin{4};
    PHI = varargin{5};
    THETA = varargin{6};
    PSI = varargin{7};
    coTr2g = varargin{8};
elseif(nargin == 10)
% getRotation('FrOfRef',MtXYZ,PHI,THETA,PSI,coTr2g,Xp0,Yp0,Zp0,flgRtOfst)
    frOfRef = varargin{1};
    [X Y Z] = getVtFromMt(varargin{2});
    PHI = varargin{3};
    THETA = varargin{4};
    PSI = varargin{5};
    coTr2g = varargin{6};
    Xp0 = varargin{7};
    Yp0 = varargin{8};
    Zp0 = varargin{9};
    flgRtOfst = varargin{10};
elseif(nargin == 12)
% getRotation('FrOfRef',X,Y,Z,PHI,THETA,PSI,coTr2g,Xp0,Yp0,Zp0,flgRtOfst)
    frOfRef = varargin{1};    
    X = varargin{2};
    Y = varargin{3};
    Z = varargin{4};
    PHI = varargin{5};
    THETA = varargin{6};
    PSI = varargin{7};
    coTr2g = varargin{8};
    Xp0 = varargin{9};
    Yp0 = varargin{10};
    Zp0 = varargin{11};
    flgRtOfst = varargin{12};
else
    error('Invalid number of input arguments');
end
    
if(PHI<0)||(PHI>pi)
    error('PHI should be between 0 and pi');
end
if(PSI<=-pi)||(PSI>pi)
    error('PSI should be between -pi and pi');
end

if(strcmpi(frOfRef,'receiver')) % if coord in receiver system
    % convert to global system
    MtRcv = [X;Y;Z];
    MtGlb = coTr2g*MtRcv;
    Xr = MtGlb(1,:);
    Yr = MtGlb(2,:);
    Zr = MtGlb(3,:);
    if flgRtOfst    % if rotate offsets
        % Add receiver offset
        Xr = Xr+Xp0;
        Yr = Yr+Yp0;
        Zr = Zr+Zp0;
    end
else                            % if coord in global system
    if flgRtOfst
        Xr = X;
        Yr = Y;
        Zr = Z;
    else   % if not rotate offsets
        % Subtract receiver offset
        Xr = X-Xp0;
        Yr = Y-Yp0;
        Zr = Z-Zp0;
    end
end
% At this point, Xr Yr Zr should be in global coordinates
%ROTATIONS IN GLOBAL SYSTEM
% get zenith angle
[Xr Yr Zr] = rotate(Xr,Yr,Zr,PHI);
% get azimuth
[Zr Xr Yr] = rotate(Zr,Xr,Yr,THETA);

% GET NEW TRANSFORMATION MATRICES
% coTr2g = getRotation('receiver',coTr2g,PHI,THETA,0,coTr2g);
[Xt Yt Zt] = rotate(coTr2g(1,:),coTr2g(2,:),coTr2g(3,:),PHI);
[Zt Xt Yt] = rotate(Zt,Xt,Yt,THETA);
coTr2g = [Xt;Yt;Zt];
coTg2r = coTr2g';

% convert to receiver coordinate system
MtGlb = [Xr;Yr;Zr];
MtRcv = coTg2r*MtGlb;
Xr = MtRcv(1,:);
Yr = MtRcv(2,:);
Zr = MtRcv(3,:);
    
% At this point, Xr Yr Zr should be in receiver coordinates
% get tilt
% NOTE: For Tilt, the receiver co-ordinates do not change
% For Tilt, only the orientation of receier system changes wrt global
% system. 
% tilt the global basis vectors in -ve PSI direction.

[Zt Xt Yt] = rotate(coTr2g(:,3),coTr2g(:,1),coTr2g(:,2),-PSI);
coTr2g = [Xt Yt Zt];

% convert to global system
MtRcv = [Xr;Yr;Zr];
MtGlb = coTr2g*MtRcv;
Xr = MtGlb(1,:);
Yr = MtGlb(2,:);
Zr = MtGlb(3,:);
if ~flgRtOfst % if not rotate offsets
    % Add receiver offset
    Xr = Xr+Xp0;
    Yr = Yr+Yp0;
    Zr = Zr+Zp0;
end

if (nargout == 1)
    Vr = zeros([3 numel(X)]);
    Vr(1,:) = Xr(:);
    Vr(2,:) = Yr(:);
    Vr(3,:) = Zr(:);
    varargout{1} = squeeze(Vr);
elseif (nargout == 2)
    Vr = zeros([3 numel(X)]);
    Vr(1,:) = Xr(:);
    Vr(2,:) = Yr(:);
    Vr(3,:) = Zr(:);
    varargout{1} = squeeze(Vr);
    varargout{2} = coTr2g;
elseif (nargout == 3)
    varargout{1} = Xr;
    varargout{2} = Yr;
    varargout{3} = Zr;
elseif (nargout == 4)
    varargout{1} = Xr;
    varargout{2} = Yr;
    varargout{3} = Zr;
    varargout{4} = coTr2g;
else
    error('Invalid number of output arguments');
end
end

%% Sub function rotation
function [ar br cr] = rotate(a,b,c,t)
Rbc = sqrt((b.^2) + (c.^2));
Thbc = atan2(c,b);
ft = t + Thbc;
ar = a;
br = Rbc.*cos(ft);
cr = Rbc.*sin(ft);
end

%% Sub function
function [X Y Z] = getVtFromMt(V)
if(size(V,1) ~= 3)
    error('Input matrix must have 3 rows');
end
X = squeeze(V(1,:));
Y = squeeze(V(2,:));
Z = squeeze(V(3,:));
end































