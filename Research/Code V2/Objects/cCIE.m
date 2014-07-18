classdef cCIE
    % class to handle CIE colorimetric data
    % specify the three observer functions. Default CIE XYZ 1978
    % calculates tristimulus and color space coordinates
    properties(SetAccess = immutable)
        ob1;    % standard observer 1
        ob2;    % standard observer 2
        ob3;    % standard observer 3
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
            if isa(psd,'cCurve')
                [X,Y,Z] = obj.getTristimulusValues(psd);
                sm = X+Y+Z;
                c1 = X./sm;
                c2 = Y./sm;
                c3 = Z./sm;
            else
                error('''psd'' argument must be of type cCurve');
            end
            
            if nargout == 0
                imshow('figCIEXYZ.png');
                axis tight;
                for iP = 1:numel(psd)
%                     x = 0.195 + c1(iP)*(0.835-0.195)/0.8;
%                     y = 0.165 + c2(iP)*(0.915-0.165)/0.9;
                    x = 32 + c1(iP)*(407-32)/0.8;
                    y = 435 + c2(iP)*(12-435)/0.9;
                    ah = annotation('ellipse','EdgeColor','k','LineWidth',2);
                    set(ah,'parent',gca);
                    set(ah,'position',[x y 10 10]);
                end
            end
        end % end getCoordinates
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
































