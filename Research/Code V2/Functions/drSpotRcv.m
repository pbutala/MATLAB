%% Comment
% drSpotRcv(Xs,Ys,as,Xrx,Yrx,arx,alrx,varargin)
% Draws the received spots on the sensor.
%
% H: handle of figure to draw the spots on.
%
% Xs Ys Zs as: Vector or Matrix specifying the centre and side length for
% each spot.
%
% Xrx Yrx arx alrx: Scalars specifying the center and side length of the
% sensor and pixel side length
%
% Specify optional input arguments in any order
%
% 1. ColorSpot: Generates the color of the spot based on the current
% colormap.
%
%% Function
function drSpotRcv(Xs1,Ys1,Xs2,Ys2,Xs3,Ys3,Xs4,Ys4,...
    pxX1,pxY1,pxX2,pxY2,pxX3,pxY3,pxX4,pxY4,varargin)
% default initialization for varargin params
Clr = 0.5;
shAll = 0;
clrPx = [1 1 1];
alPx = 0.25;
% varargin check
nVArg = nargin-16;
if (rem(nVArg,2)~= 0)
    error('Check input arguments');
end
ArgName = 1;
ArgParam = ArgName + 1;
while(ArgName < nVArg)
    switch lower(varargin{ArgName})
        case 'colorspot'
            Clr = varargin{ArgParam};
        case 'showspot'
            switch lower(varargin{ArgParam})
                case 'all'
                    shAll = true;
                case 'onsensor'
                    shAll = false;
                otherwise
                    error('''Show'' parameter should be either ''All'' or ''OnSensor''');
            end
        case 'colorpixel'
            clrPx = varargin{ArgParam};
        case 'alphapixel'
            alPx = varargin{ArgParam};
        otherwise
            error('unknown parameter %s specified',varargin{ArgName});
    end
    ArgName = ArgName + 2;
    ArgParam = ArgName + 1;
end

if(~isequal(size(Xs1),size(Ys1),size(Xs2),size(Ys2),...
        size(Xs3),size(Ys3),size(Xs4),size(Ys4)))
    error('All Xs Ys dimensions must agree');
end

% if(~isscalar(Xrx))
%     error('Xrx Yrx arx alrx must be scalar.');
% end
% 
% if(~isequal(size(Xrx),size(Yrx),size(arx),size(alrx)))
%     error('Xrx Yrx arx alrx dimensions must agree');
% end

if(numel(Clr)==1)
    tClr = repmat(Clr,size(Xs1));
elseif(isequal(size(Clr),size(Xs1)))
    tClr = Clr;
else
    error('ColorSpot must be a scalar or a matrix with same dimensions as Xs and Ys');
end

if size(clrPx,1) == 1
    clrPx = repmat(clrPx,numel(pxX1),1);
end

if isscalar(alPx)
    alPx = repmat(alPx,size(pxX1));
end
%% Logic
Ns = numel(Xs1);
Xs = [pxX1(:);pxX2(:);pxX3(:);pxX4(:)];     
gdXmin = min(Xs);
gdXmax = max(Xs);
arx = gdXmax - gdXmin;
% Ys = [pxY1(:);pxY2(:);pxY3(:);pxY4(:)];
% gdYmin = min(Ys);
% gdYmax = max(Ys);
% dY = gdYmax - gdYmin;

for idPx = 1:numel(pxX1)
    vtx1 = [pxX1(idPx),pxY1(idPx);...
        pxX2(idPx),pxY2(idPx);...
        pxX3(idPx),pxY3(idPx);...
        pxX4(idPx),pxY4(idPx)];
     patch('Vertices',vtx1,'Faces',[1 2 3 4 1],'EdgeColor',[0 0 0],'FaceColor',clrPx(idPx,:),'FaceAlpha',alPx(idPx));
end

drx = arx/2;
is = 1;
while(is <= Ns)
    TRX = Xs1(is);
    TLX = Xs2(is);
    BLX = Xs3(is);
    BRX = Xs4(is);
    
    TRY = Ys1(is);
    TLY = Ys2(is);
    BLY = Ys3(is);
    BRY = Ys4(is);
    
    if(~shAll)
        TRX(TRX>drx) = drx;
        TRX(TRX<-drx) = -drx;
        TRY(TRY>drx) = drx;
        TRY(TRY<-drx) = -drx;
        
        TLX(TLX>drx) = drx;
        TLX(TLX<-drx) = -drx;
        TLY(TLY>drx) = drx;
        TLY(TLY<-drx) = -drx;
        
        BLX(BLX>drx) = drx;
        BLX(BLX<-drx) = -drx;
        BLY(BLY>drx) = drx;
        BLY(BLY<-drx) = -drx;
        
        BRX(BRX>drx) = drx;
        BRX(BRX<-drx) = -drx;
        BRY(BRY>drx) = drx;
        BRY(BRY<-drx) = -drx;
        
        set(gca,'XLimMode','manual','XLim',[-arx/2 arx/2]);
        set(gca,'YLimMode','manual','YLim',[-arx/2 arx/2]);
    end
    
    % Draw spots
    TRs = [TRX TRY];
    TLs = [TLX TLY];
    BLs = [BLX BLY];
    BRs = [BRX BRY];
    vtx = [TRs;TLs;BLs;BRs];
    fcs = [1 2 3 4];
    patch('Vertices',vtx,'Faces',fcs,'FaceVertexCData',tClr(is),'FaceColor','flat');
      
    is = is+1;
end

% Draw the sensor grid
% for xp=-drx:alrx:drx
%     vtx1 = [drx xp;-drx xp];
%     patch('Vertices',vtx1,'Faces',[1 2],'EdgeColor',[0 0 0]);
% end
% for yp=-drx:alrx:drx
%     vtx1 = [yp drx;yp -drx];
%     patch('Vertices',vtx1,'Faces',[1 2],'EdgeColor',[0 0 0]);
% end
% for idPx = 1:numel(pxX1)
%     vtx1 = [pxX1(idPx),pxY1(idPx);...
%         pxX2(idPx),pxY2(idPx);...
%         pxX3(idPx),pxY3(idPx);...
%         pxX4(idPx),pxY4(idPx)];
%      patch('Vertices',vtx1,'Faces',[1 2 3 4 1],'EdgeColor',[0 0 0],'FaceColor',[1 1 1],'FaceAlpha',0.25);
% end

% set(0,'CurrentFigure',CF);

%% Space
































