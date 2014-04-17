classdef cImagingReceiver
    properties
        location; % location of receiver in global co-ordinates
        orientation = cOrientation(0,0,0); % orientation of receiver
    end
    properties(SetAccess = private)
        optics; % imaging
        sensor; % imaging sensor
        tia;    % transimpedance amplifier
        % combining
    end
    properties(Dependent = true, SetAccess = private)
        rxCount;    % total number of receivers of this type
        rxFOV;      % field of view of receiver
    end
    
    methods
        % Constructor
        function obj = cImagingReceiver(Xr,Yr,Zr,Xs,Ys,f,fN)
            % cImagingReceiver(Xr,Yr,Zr,Xs,Ys,f,fN): Class constructor
            % Xr,Yr,Zr: X,Y,Z coord of receiver locations. Default (0,0,0)
            % Xs,Ys: X,Y coord of center of pixel in local coord
            % f,fN: focal length, aperture f/#
            obj.location = cLocation(Xr,Yr,Zr);
            obj.optics = cImagingOptics(f,fN);
            obj.sensor = cSensor(Xs,Ys,repmat(-f,size(Xs)));
            obj.tia = cTIA(5e-12*ones(size(Xs)));
        end
    end
    
    methods
        function [spX1 spY1 spX2 spY2 spX3 spY3 spX4 spY4] =...
                getImage(obj,rxLoc,txLoc,txSz,txOri,Hax)
            rxOri = obj.orientation;
            % Calculate Received Image
            [spX1 spY1 spX2 spY2 spX3 spY3 spX4 spY4] =...
                getRcvImg(txLoc.X,txLoc.Y,txLoc.Z,...
                rxLoc.X,rxLoc.Y,rxLoc.Z,...
                txSz.L,obj.optics.f,...
                'ZenithTx',txOri.Z(:),...
                'AzimuthTx',txOri.A(:),...
                'TiltTx',txOri.T(:),...
                'ZenithRx',rxOri.Z(:),...
                'AzimuthRx',rxOri.A(:),...
                'TiltRx',rxOri.T(:));
            % if no output arguments, draw the image on current axis
            if nargout == 0
                    pxClr = zeros(numel(obj.sensor.pxR),3);
                    for idPx = 1:numel(obj.sensor.pxR)
                        pxClr(idPx,:) = obj.sensor.color(obj.sensor.type(idPx),:);
                    end
                for idRx = 1:numel(rxLoc.X)
                    if exist('Hax','var')
%                         axes(Hax(idRx));
                        set(0,'CurrentFigure',get(Hax(idRx),'Parent'));
                    else
                        figure;
                    end
                    drSpotRcv(spX1(idRx,:),spY1(idRx,:),...
                        spX2(idRx,:),spY2(idRx,:),...
                        spX3(idRx,:),spY3(idRx,:),...
                        spX4(idRx,:),spY4(idRx,:),...
                        obj.sensor.pxR,obj.sensor.pxT, obj.sensor.pxL,obj.sensor.pxT,...
                        obj.sensor.pxL,obj.sensor.pxB, obj.sensor.pxR,obj.sensor.pxB,...
                        'Showspot','OnSensor','ColorPixel',pxClr,'AlphaPixel',0.2);
                    
                    % cb = colorbar;
                    % set(cb,'ylim',[cmin cmax]/1e6, 'ytick',cmin/1e6:(cmax-cmin)/1e7:cmax/1e6);
                    
                    xlabel('Sensor Length (mm)');
                    ylabel('Sensor Width (mm)');
                    tStr = sprintf('Image formed on sensor by transmitters\nReceiver location = [%0.2f %0.2f %0.2f]',rxLoc.X,rxLoc.Y,rxLoc.Z);
                    title(tStr);
                    rxa = obj.sensor.edge.L;
                    ctTickRes = 1e-3;
                    sTicks = -(rxa/2)-rem(-(rxa/2),ctTickRes):ctTickRes:(rxa/2);
%                     sTicksLbl = sprintf('%0.3f|',sTicks(:));
                    set(gca,'XTickMode','Manual','XTickLabelMode','Manual');
                    set(gca,'YTickMode','Manual','YTickLabelMode','Manual');
                    set(gca,'XTickLabel',[]);
                    set(gca,'YTickLabel',[]);
                    set(gca,'XTick',sTicks);
                    set(gca,'YTick',sTicks);
%                     set(gca,'XTick',sTicks,'XTickLabel',sTicksLbl);
%                     set(gca,'YTick',sTicks,'YTickLabel',sTicksLbl);
%                     set(gca,'XTick',sTicks);
%                     set(gca,'XTickLabel',{'one';'two'});
%                     set(gca,'YTick',sTicks);
%                     set(gca,'YTickLabel',{'one';'two'});
                    axis equal;
                    axis tight;
                end
            end
        end
    end
    methods(Abstract)
        [Ipx Iambpx Isig Iamb frSpOnPx frPxOnSp] = ...
                getSignal(obj,rxPSD,ambPSD,txLoc,txSz,txOri);
    end
    
    %% Get/Set accessors
    methods
        function n = get.rxCount(obj)
            n = obj.location.count;
        end
        
        function fov = get.rxFOV(obj)
            % gets the field of view half angle of the receiver. depends on
            % the sensor size and the fov of the optics used
            diag = sqrt((obj.sensor.edge.L^2)+(obj.sensor.edge.W^2));
            fov = atan(diag/(2*obj.optics.f));
            if fov > obj.optics.fov
                fov = obj.optics.fov;
            end
        end
    end
end