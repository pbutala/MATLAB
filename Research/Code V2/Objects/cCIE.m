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
        end
        
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
        end
        
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
        end
    end
end
































