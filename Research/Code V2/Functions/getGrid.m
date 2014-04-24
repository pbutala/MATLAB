function [Xs,Ys,Zs] = getGrid(varargin)
% Computes and returns grid coordinates
% [Xs Ys Zs] = getGrid(X,Y,Z)
% -INPUT-
% X,Y,Z: number of elements
% -OUTPUT-
% Xs Ys Zs: centered grid about 0.
%
% [Xs Ys Zs] = getGrid(X,Y,Z,Px,Py,Pz)
% -INPUT-
% X,Y,Z: number of elements
% Px,Py,Pz: element pitch
% -OUTPUT-
% Xs Ys Zs: centered grid about 0.
%
% [Xs Ys Zs] = getGrid(X,Y,Z,Px,Py,Pz,'Fill')
% -INPUT-
% X,Y,Z: cube dimensions
% Px,Py,Pz: element pitch
% -OUTPUT-
% Xs Ys Zs: centered grid about 0.
X = varargin{1};Y = varargin{2};Z = varargin{3};
Px = 1;Py = 1;Pz = 1;
if nargin > 3
    Px = varargin{4};Py = varargin{5};Pz = varargin{6};
    if strcmpi(varargin{end},'Fill')
        if (Px == 0) || (Py == 0) || (Pz == 0)
            error('Pitch cannot be ''0'' for ''Fill''')
        end
        X = floor(X/Px)+1;Y = floor(Y/Py)+1;Z = floor(Z/Pz)+1;
    end
end
Xa = getArray(X,Px);
Ya = getArray(Y,Py);
Za = getArray(Z,Pz);

[Xs,Ys,Zs] = meshgrid(Xa,Ya,Za);
end