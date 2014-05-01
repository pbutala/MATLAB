function drCrosshair(X,Y,fMrk,type,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isequal(size(X),size(Y))
    error('X and Y must have same lengths');
end
if ~exist('lSty','var')
    lSty = 'k-';
end
if ~exist('fMrk','var')
    fMrk = true;
end
if ~isa(fMrk,'logical')
    error('fMRK, must be a logical- ''TRUE'' or ''FALSE''');
end
if ~(strcmpi(type,'full')||strcmpi(type,'half'))
    error('Type must be either ''FULL'' or ''HALF''');
end
NP = numel(X);
% if numel(lSty)==1
%     lSty = repmat(lSty,NP,1);
% end
% if size(lSty,1)~=NP
%     error('Number of linestyles must be equal to either 1 or number of points');
% end

if numel(fMrk)==1
    fMrk = repmat(fMrk,NP,1);
end
if numel(fMrk)~=NP
    error('Number of marker flags must be equal to either 1 or number of points');
end

AX = gca;
XLIM = get(AX,'XLim');
YLIM = get(AX,'YLim');
S = get(0,'DefaultLineMarkerSize')*4;
C = [0 0 0];
hold all;
for i=1:NP
    % Get point, linestyle and markerflag
    Px = X(i);Py = Y(i);
    LS = lSty(i,:);
    FMRK = fMrk(i);
    % Draw line parallel to x-axis through (X,Y)
    if strcmpi(type,'half')
        plot([XLIM(1) Px],[Py Py],varargin{:});
    else
        plot(XLIM,[Py Py],varargin{:});
    end
    % Draw line parallel to y-axis through (X,Y)
    if strcmpi(type,'half')
        plot([Px Px],[YLIM(1) Py],varargin{:});
    else
        plot([Px Px],YLIM,varargin{:});
    end
    % Draw marker if fMrk == true
    if(FMRK)
        scatter(Px,Py,S,C);
    end
end
end





























