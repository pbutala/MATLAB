% cRoom: class to handle a Room
classdef cRoom < handle
    properties
        luminaire;  % collection of all luminaire types in the room
        % additional isotropic ambient light psd
        % rcvLocation;
        % walls
        % reflectivities
        % etc ROOM specific params
    end
    properties(Access = private)
        dimension;  % dimension of the room
    end
    properties(Dependent = true, SetAccess = private)
        L;              % gets the length of the room
        W;              % gets the width of the room
        H;              % gets the height of the room
        lmCount;        % gets number of luminaires in the room
        maxLmCount;     % gets maximum number of luminaires of a luminaire type
        chCount;        % gets maximum number of channels for each luminaire
        maxChCount;     % gets maximum number of channels that a luminaire may have  
        lmLocation;     % gets locations of all transmitters in the room
        lmDimension;    % gets dimension of all transmitters in the room
    end
    
    %% Class methods
    methods
        function obj = cRoom(o1,o2,o3)
            % cRoom(L,W,H): class constructor
            % cRoom(cSize): class constructor
            % -INPUT- 
            % L,W,H: length,width,hright of room. (default 7,4,4)
            if nargin == 1
                obj.dimension = o1;
            elseif nargin > 1
                obj.dimension = cSize(o1,o2,o3);
            else
                obj.dimension = cSize(7,4,4);
            end
        end
        
        function H = getFreeSpaceGain(obj,rxLoc,rxOri,rxFOV)
            % H = getFreeSpaceGain(rxLoc)
            % Calculates and returns the free space gain 
            % -INPUT-
            % rxLoc: receiver locations to calculate gain to
            % -OUTPUT- 
            % H: gain matrix (rx,lmtyp,lm,ch)
            
            % if takes long time, may want to pre-compute H and store in
            % memory
            nRx = numel(rxLoc.X);
            H = zeros(nRx,numel(obj.luminaire),obj.maxLmCount,obj.maxChCount);
            nL = numel(obj.luminaire);
            for iL = 1:nL
                nTx = obj.luminaire(iL).lmCount;
                nC = obj.luminaire(iL).chCount;
                H(:,iL,1:nTx,1:nC) = ...
                    obj.getFreeSpaceGainLuminaireChannels(obj.luminaire(iL),rxLoc,rxOri,rxFOV);
            end
        end
        
        function H = getFreeSpaceGainLuminaireChannels(obj,lmn,rxLoc,rxOri,rxFOV)
            % H = getFreeSpaceGainLuminaireChannels(lmn,rxLoc)
            % Calculates and returns the free space gain
            % -INPUT-
            % lmn: luminaire to calculate free space gain for
            % rxLoc: receiver locations to calculate gain to
            % -OUTPUT- 
            % H: gain matrix(rx,tx,ch)
            
            % TODO: account for reflections (rx,lm,ch)
        
            % if takes long time, may want to pre-compute H and store in
            % memory
            nRx = numel(rxLoc.X);
            nTx = lmn.lmCount;
            nC = lmn.chCount;
            H = zeros(nRx,nTx,nC);
            for iC = 1:nC
                % if not, create and store
                [PHI,PSI,D] = getLOSParams(lmn.location.X(:),lmn.location.Y(:),...
                    lmn.location.Z(:),rxLoc.X(:),rxLoc.Y(:),rxLoc.Z(:),...
                    'ZenithTx',lmn.orientation.Z(:),...
                    'AzimuthTx',lmn.orientation.A(:),...
                    'TiltTx',lmn.orientation.T(:),...
                    'ZenithRx',rxOri.Z(:),...
                    'AzimuthRx',rxOri.A(:),...
                    'TiltRx',rxOri.T(:));
                h = ((lmn.order + 1)*(cos(PHI).^lmn.order).*cos(PSI))./(2*pi*(D.*D));
                h(PSI>rxFOV) = 0; % set gain = 0 for locations outside FOV
                H(:,:,iC) = h;
                % for LOS each color has same channel gain. 
                % when colored reflections are considered, different color
                % channels will have different gains.
                % 'obj' will have reflectivity information.
            end
        end
        
        function IrdClr = getIrradianceLuminaireChannels(obj,lmn,rxLoc,rxOri,rxFOV)
            % IrdClr = getIrradianceLuminaireChannels(lmn,rxLoc)
            % gets irradiance for each color channel from each location of the
            % luminaire type
            % -INPUT-
            % lmn: luminaire to calculate free space gain for
            % rxLoc: receiver locations to calculate gain to
            % -OUTPUT- 
            % IrdClr: Irradiance(rx,tx,ch)
            
            H = obj.getFreeSpaceGainLuminaireChannels(lmn,rxLoc,rxOri,rxFOV);
            IrdClr = zeros(size(H));
            fClr = lmn.rdFluxClr;
            nC = numel(fClr);
            for iC=1:nC
                IrdClr(:,:,iC) = H(:,:,iC)*fClr(iC);
            end
        end
        
        function IlmClr = getIlluminanceLuminaireChannels(obj,lmn,rxLoc,rxOri,rxFOV)
            % IlmClr = getIlluminanceLuminaireChannels(lmn,rxLoc)
            % gets illuminance for each color channel from each location of the
            % luminaire type
            % -INPUT-
            % lmn: luminaire to calculate free space gain for
            % rxLoc: receiver locations to calculate gain to
            % -OUTPUT- 
            % IlmClr: Illuminance(rx,tx,ch)
            H = obj.getFreeSpaceGainLuminaireChannels(lmn,rxLoc,rxOri,rxFOV);
            IlmClr = zeros(size(H));
            fClr = lmn.lmFluxClr;
            nC = numel(fClr);
            for iC=1:nC
                IlmClr(:,:,iC) = H(:,:,iC)*fClr(iC);
            end
        end
        
        function Irad = getIrradiance(obj,rxLoc,rxOri)
            % Irad = getIrradiance(rxLoc)
            % gets the total irradiance
            % -INPUT-
            % rxLoc: receiver locations
            % -OUTPUT- 
            % Irad: Irradiance(rx,tx,ch)
            Irad = zeros(numel(rxLoc.X),1);
            nL = numel(obj.luminaire);
            for iL = 1:nL
                IrdC = obj.getIrradianceLuminaireChannels(obj.luminaire(iL),rxLoc,rxOri,pi/2);
                Irad = Irad + sum(sum(IrdC,3),2);
            end
        end
        
        function Ilm = getIlluminance(obj,rxLoc,rxOri)
            % Ilm = getIlluminance(rxLoc)
            % gets the total illuminance
            % -INPUT-
            % rxLoc: receiver locations
            % -OUTPUT- 
            % Ilm: Illuminance(rx,tx,ch)
            Ilm = zeros(numel(rxLoc.X),1);
            nL = numel(obj.luminaire);
            for iL = 1:nL
                IlmC = obj.getIlluminanceLuminaireChannels(obj.luminaire(iL),rxLoc,rxOri,pi/2);
                Ilm = Ilm + sum(sum(IlmC,3),2);
            end
        end
        
        function setIlluminance(obj,rxLoc,vIll)
            % the function calculates the scale to multiply the luminaire
            % output powers with such that the rxLoc specified will have
            % vIll illumination
            % rxLoc must point to one location only.
            Ill = obj.getIlluminance(rxLoc,cOrientation(0,0,0));
            scale = vIll/Ill;
            nL = numel(obj.luminaire);
            for iL = 1:nL
                obj.luminaire(iL).scaleOutputFlux(scale);
            end
        end
        
    end
    
    %% drawing methods
    methods
        function Ilm = drawIlluminance(obj,plLoc,ilTh,Hax)
            % H = drawIlluminance(obj,plLoc)
            % Draws illuminance on the current axis
            % -INPUT-
            % plLoc: plane locations
            % -OUTPUT-
            % H: handle to the surface
            if nargin == 1 % plLoc not provided
                [plX plY plZ] = getGrid(obj.L,obj.W,1,0.2,0.2,2,'Fill');
                plX = plX + obj.L/2;
                plY = plY + obj.W/2;
                plZ = plZ + 1;
                plLoc = cLocation(plX,plY,plZ);
            end
            Ilm = obj.getIlluminance(plLoc,cOrientation(0,0,0));
            if exist('Hax','var')
                H = Hax(1);
            else
                figure;
                H = gca;
            end
            set(0,'CurrentFigure',get(H,'Parent'));
            surf(plLoc.X,plLoc.Y,reshape(Ilm,size(plLoc.X)),'FaceColor','interp');
            minY = min(Ilm(:));
            maxY = max(Ilm(:));
            if maxY == minY,
                maxY = maxY + 1;
            end
            axis([0 obj.L 0 obj.W minY maxY]);
            grid on;
            colorbar;
            view(3);
            tStr = sprintf('Max: %0.2f lx',max(Ilm(:)));
            title(tStr);
            xlabel('X');
            ylabel('Y');
            zStr = sprintf('Illuminance (lx)');
            zlabel(zStr);
            if exist('ilTh','var')
                IlmCvg = zeros(size(Ilm));
                IlmCvg(Ilm >= ilTh) = 1;
                ilCvgPerc = (sum(IlmCvg(:))*100)/numel(IlmCvg);
                if exist('Hax','var')
                    H = Hax(2);
                else
                    figure;
                    H = gca;
                end
                set(0,'CurrentFigure',get(H,'Parent'));
                image(plLoc.X(1,:),plLoc.Y(:,1),reshape(100*IlmCvg,size(plLoc.X)));
%                 axis([0 obj.L 0 obj.W min(Ilm(:)) max(Ilm(:))]);
%                 grid on;
                tStr = sprintf('Illuminance Coverage\n%0.2f%%',ilCvgPerc);
                title(tStr);
                xlabel('X');
                ylabel('Y');
                zStr = sprintf('Illuminance Coverage (%%)');
                zlabel(zStr);
            end
        end
        
        function Irad = drawIrradiance(obj,plLoc)
            % H = drawIrradiance(obj,plLoc)
            % Draws irradiance on the current axis
            % -INPUT-
            % plLoc: plane locations
            % -OUTPUT-
            % H: handle to the surface
            if nargin == 1 % rxLoc not provided
                [plX plY plZ] = getGrid(obj.L,obj.W,1,0.2,0.2,2,'Fill');
                plX = plX + obj.L/2;
                plY = plY + obj.W/2;
                plZ = plZ + 1;
                plLoc = cLocation(plX,plY,plZ);
            end
            Irad = obj.getIrradiance(plLoc,cOrientation(0,0,0));
            surf(plLoc.X,plLoc.Y,reshape(Irad,size(plLoc.X)),'FaceColor','interp');
            minY = min(Irad(:));
            maxY = max(Irad(:));
            if maxY == minY,
                maxY = maxY + 1;
            end
            axis([0 obj.L 0 obj.W minY maxY]);
            grid on;
            colorbar;
            view(3);
            tStr = sprintf('Irradiance\nMax: %0.2f W.m^{-2}',max(Irad(:)));
            title(tStr);
            xlabel('X');
            ylabel('Y');
        end
        
        function drawSetup(obj,rxLoc,rxOri,rxFOV)
            % H = drawRoomSetup(rxLoc)
            % draws the room setup on current axis
            txLoc = obj.lmLocation;
            txOri = obj.luminaire.orientation;
            drSetupImg(obj.L,obj.W,obj.H,...
                txLoc.X,txLoc.Y,txLoc.Z,...
                rxLoc.X,rxLoc.Y,rxLoc.Z,'LightRays','FOV',...
                'FOVhalf',rxFOV,...
                'ZenithTx',txOri.Z(:),...
                'AzimuthTx',txOri.A(:),...
                'TiltTx',txOri.T(:),...
                'ZenithRx',rxOri.Z(:),...
                'AzimuthRx',rxOri.A(:),...
                'TiltRx',rxOri.T(:));
%             ,...
%                 'FOVhalf',Vho,...
%                 'LightRays','FOV',...
%                 'ZenithTx',txZeta,'AzimuthTx',txAlpha,'TiltTx',txTau,...
%                 'ZenithRx',rxZeta,'AzimuthRx',rxAlpha,'TiltRx',rxTau
            axis([0 obj.L 0 obj.W 0 obj.H]);
            
            set(gca,'XTick',0:1:obj.L);
            set(gca,'YTick',0:1:obj.W);
            set(gca,'ZTick',0:1:obj.H);
            
            grid on;
            view(3);
            xlabel('X: Length (m)');
            ylabel('Y: Width (m)');
            zlabel('Z: Height (m)');
            title('System Setup');
        end
    end
    
    %% Property Get/Set accessors
    methods
        function L = get.L(obj)
            % gets the length of the room
            L = obj.dimension.L;
        end
        
        function W = get.W(obj)
            % gets the width of the room
            W = obj.dimension.W;
        end
        
        function H = get.H(obj)
            % gets the height of the room
            H = obj.dimension.H;
        end
        
        function n = get.lmCount(obj)
            % gets number of luminaires in the room
            nL = numel(obj.luminaire);
            n = zeros(size(obj.luminaire));
            for iL = 1:nL
                n(iL) = obj.luminaire(iL).lmCount;
            end
        end
        
        function n = get.maxLmCount(obj)
            % gets maximum number of luminaires of a luminaire type
            n = max(obj.lmCount(:));
        end
        
        function n = get.chCount(obj)
            % gets maximum number of channels for each luminaire
            nL = numel(obj.luminaire);
            n = zeros(size(obj.luminaire));
            for iL = 1:nL
                n(iL) = obj.luminaire(iL).chCount;
            end
        end
         
        function n = get.maxChCount(obj)
            % gets maximum number of channels that a luminaire may have
            n = max(obj.chCount(:));
        end
        
        function Loc = get.lmLocation(obj)
            % gets locations of all transmitters in the room
            nL = numel(obj.luminaire);
            Loc = obj.luminaire(1).location;
            for iL = 2:nL
                Loc = Loc + obj.luminaire(iL).location;
            end
        end
        
        function Dim = get.lmDimension(obj)
            % gets dimension of all transmitters in the room
            nL = numel(obj.luminaire);
            Dim = obj.luminaire(1).dimension;
            for iL = 2:nL
                Dim = Dim + obj.luminaire(iL).dimension;
            end
        end
    end
end






























