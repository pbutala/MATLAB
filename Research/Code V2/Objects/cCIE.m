classdef cCIE < handle
    % class to handle CIE colorimetric data
    % specify the three observer functions. Default CIE XYZ 1978
    % calculates tristimulus and color space coordinates
    properties(SetAccess = immutable)
        ob1;    % standard observer 1
        ob2;    % standard observer 2
        ob3;    % standard observer 3
    end
    
    properties (SetAccess = private)
        gmtX;
        gmtY;
        gmtL;
    end
    
    methods
        function x = get.gmtX(obj)
            if isempty(obj.gmtX)
                [obj.gmtX,obj.gmtY] = obj.getGamut();
            end
            x = obj.gmtX;
        end
        
        function y = get.gmtY(obj)
            if isempty(obj.gmtY)
                [obj.gmtX,obj.gmtY] = obj.getGamut();
            end
            y = obj.gmtY;
        end
        
        function l = get.gmtL(obj)
            if isempty(obj.gmtL)
                [obj.gmtX,obj.gmtY,obj.gmtL] = obj.getGamut();
            end
            l = obj.gmtL;
        end
    end
    
    methods
        function obj = cCIE(file)
            % Class constructor
            if nargin == 0
                file = 'CIE1978.csv';
            end
            fl = load(file);
            L = fl(:,1);
            X = fl(:,2);
            Y = fl(:,3);
            Z = fl(:,4);
            Lmin = L(1);
            Lmax = L(end);
            dL = L(2)-L(1);
            obj.ob1 = cCurve(Lmin,dL,Lmax,X);
            obj.ob2 = cCurve(Lmin,dL,Lmax,Y);
            obj.ob3 = cCurve(Lmin,dL,Lmax,Z);
            obj.gmtX = [];
            obj.gmtY = [];
        end % end constructor
        
        function [t1,t2,t3] = getTristimulusValues(obj,psd)
            % [t1 t2 t3] = getTristimulusValues(psd)
            % computes and returns the tristimulus values for the psd
            if isa(psd,'cCurve')
                t1 = zeros(size(psd));
                t2 = zeros(size(psd));
                t3 = zeros(size(psd));
                for iP = 1:numel(psd)
                    v = psd(iP).*obj.ob1;
                    t1(iP) = v.getIntegral(380,780);
                    v = psd(iP).*obj.ob2;
                    t2(iP) = v.getIntegral(380,780);
                    v = psd(iP).*obj.ob3;
                    t3(iP) = v.getIntegral(380,780);
                end
            else
                error('''psd'' argument must be of type cCurve');
            end
        end % end getTristimulusValues
        
        function [c1,c2,c3] = getCoordinates(obj,psd)
            % [c1 c2 c3] = getCoordinates(psd)
            % computes and returns the CIE xyz valuescolor space coordinates
            % for the psd
            if isa(psd,'cPSD')
                [X,Y,Z] = obj.getTristimulusValues(psd.npsd);
            elseif isa(psd,'cCurve')
                [X,Y,Z] = obj.getTristimulusValues(psd);
            elseif isscalar(psd) %  assume is wavelength
                psd = floor(psd);
                if((psd<380) || (psd>780))
                    error('Wavelength must be between 380nm and 780nm');
                end
                Ys = zeros(1,401);
                Ys(psd-379) = 1;
                c = cCurve(380,1,780,Ys);
                [X,Y,Z] = obj.getTristimulusValues(c);
            else
                error('''psd'' argument must be of type cCurve');
            end
            sm = X+Y+Z;
            c1 = X./sm;
            c2 = Y./sm;
            c3 = Z./sm;
                
            if nargout == 0
                error('no output arguments specified');
            end
        end % end getCoordinates
        
        function [x,y,l] = getGamut(obj,dL,LMIN,LMAX)
            % function getGamut()
            % shows color gamut on plot
            if ~exist('dL','var')
                dL = 5;
            end
            if ~exist('LMIN','var')
                LMIN = 380;
            end
            if ~exist('LMAX','var')
                LMAX = 780;
            end
            
            l = LMIN:dL:LMAX;
            lenL = length(l);
            x = zeros(1,lenL);
            y = zeros(1,lenL);
            
            for il = 1:lenL
                Ys = zeros(1,lenL);
                Ys(il) = 1;
                psd = cCurve(LMIN,dL,LMAX,Ys);
                [x(il),y(il),~] = obj.getCoordinates(psd);
            end
        end
        
        function hax = showGamut(obj,fLabel)
            hax = gca;
            plot(gca,[obj.gmtX obj.gmtX(1)],[obj.gmtY obj.gmtY(1)],'k','LineWidth',2);
            axis([0 0.8 0 0.9]);
            axis square;
            grid on;
            if ~exist('fLabel','var')
                fLabel = true;
            end
            if fLabel
                LDATXT = [380;470;490;500;510;520;540;560;580;600;780];
                for iL = 1:length(obj.gmtL)
                    if ~isempty(find(LDATXT == obj.gmtL(iL),1))
                        str = sprintf('  %d',obj.gmtL(iL));
                        text(obj.gmtX(iL),obj.gmtY(iL),str);
                    end
                end
            end
            xlabel('x'); ylabel('y');
        end
        
        function hax = showGamutPlanckian(obj,CCT)
            hax = obj.showGamut();
            hold on;
            if ~exist('CCT','var')
                CCT = [1667 1750:250:25000];
            end
            [x,y] = planckXY(CCT);   % Get x,y from CCT
            plot(hax,x,y,'k','LineWidth',2);
        end
        
    end % end methods
    
    methods(Static)
        function varargout = xyz2rgb(varargin)
            % RGB = xyz2rgb(XYZ)
            % [R G B] = xyz2rgb([X Y Z])
            switch nargin
                case 1
                    tmp = varargin{1};
                    XYZ = tmp(:);
                case 3
                    XYZ = [varargin{1};varargin{2};varargin{3}];
                otherwise
                    error('Incorrect input arguments. Call xyz2rgb(XYZ) OR xyz2rgb(X,Y,Z)');
            end
            Trgb2xyz =  (1/0.17697)*[0.49 0.31 0.20; 0.17697 0.81240 0.01063; 0.00 0.01 0.99];
            RGB = Trgb2xyz\XYZ;
            switch nargout
                case 1
                    varargout{1} = RGB;
                case 3
                    varargout{1} = RGB(1);
                    varargout{2} = RGB(2);
                    varargout{3} = RGB(3);
                otherwise
                    error('Incorrect output arguments.');
            end
        end
        
        function varargout = rgb2xyz(varargin)
            % XYZ = rgb2xyz(RGB)
            % [X Y Z] = rgb2xyz([R G B])
            switch nargin
                case 1
                    tmp = varargin{1};
                    RGB = tmp(:);
                case 3
                    RGB = [varargin{1};varargin{2};varargin{3}];
                otherwise
                    error('Incorrect input arguments. Call rgb2xyz(RGB) OR rgb2xyz(R,G,B)');
            end
            Trgb2xyz =  (1/0.17697)*[0.49 0.31 0.20; 0.17697 0.81240 0.01063; 0.00 0.01 0.99];
            XYZ = Trgb2xyz*RGB;
            switch nargout
                case 1
                    varargout{1} = XYZ;
                case 3
                    varargout{1} = XYZ(1);
                    varargout{2} = XYZ(2);
                    varargout{3} = XYZ(3);
                otherwise
                    error('Incorrect output arguments.');
            end
        end
    end % end methods(Static)
end % end classdef
































