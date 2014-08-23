classdef cCurve < handle
    properties(SetAccess = private)
        Xmin;
        dX;
        Xmax;
        Y;
    end
    
    properties
        Ylow = 0;
        Yhigh = 0;
    end
    
    properties (Dependent = true, SetAccess = private)
        X;
        Ymin;
        Ymax;
    end
    
    methods
        % Constructor
        function obj = cCurve(Xmin,dX,Xmax,Y)
            if nargin ~= 0
                obj.initCurve(Xmin,dX,Xmax,Y);
            end
        end
        
        % linear interpolation
        function ys = getY(obj,xs)
            Xs = obj.X;
            Ys = obj.Y;
            ys = zeros(size(xs));
            for i = 1:numel(xs)
                x = xs(i);
                u = find(Xs>=x,1,'first');
                
                if isempty(u) % if u is [], u==1 check throws error
                    % All Xs < x
                    ys(i) = obj.Yhigh;
                elseif u == 1
                    % All Xs > x
                    if(Xs(1) == x)
                        ys(i) = Ys(1);
                    else
                        ys(i) = obj.Ylow;
                    end
                else
                    l = u-1;
                    m = (Ys(u)-Ys(l))/(Xs(u)-Xs(l));
                    ys(i) = Ys(l) + m*(x-Xs(l));
                end
                
            end
        end
        
        % get Integral
        function I = getIntegral(obj,xl,xu)
            I = 0;
            if nargin == 1
                xl = obj.Xmin;
                xu = obj.Xmax;
            end
            if xl < obj.Xmin
                I = I + (obj.Xmin - xl)*obj.Ymin;
                xl = obj.Xmin;
            end
            
            if xu > obj.Xmax
                I = I + (xu - obj.Xmax)*obj.Ymax;
                xu = obj.Xmax;
            end
            
            ys = obj.getY(xl:obj.dX:xu);
            I = I+sum(ys)*obj.dX;
        end
        
        % Overload '.*' operator
        function res = times(obj1,obj2)
            if ~isa(obj1,'cCurve')
                error('mtimes not defined for objects of type %s',class(obj1));
            end
            if ~isa(obj2,'cCurve')
                error('mtimes not defined for objects of type %s',class(obj2));
            end
            xmin = min(obj1.Xmin,obj2.Xmin);
            dx = min(obj1.dX,obj2.dX);
            xmax = max(obj1.Xmax,obj2.Xmax);
            Xs = xmin:dx:xmax;
            Y1s = obj1.getY(Xs);
            Y2s = obj2.getY(Xs);
            Ys = Y1s.*Y2s;
            res = cCurve(xmin,dx,xmax,Ys);
        end
        
        % Overload '*' operator
        function res = mtimes(obj1,obj2)
            if (isa(obj1,'cCurve') && isnumeric(obj2))
                crv = obj1;
                scl = obj2;
            elseif (isa(obj2,'cCurve') && isnumeric(obj1))
                crv = obj2;
                scl = obj1;
            else
                error('Input arguments must be of types numeric and cCurve');
            end
            
            res = cCurve(crv.Xmin,crv.dX,crv.Xmax,crv.Y*scl);
        end
        
        % initialize Curve
        function initCurve(obj,Xmin,dX,Xmax,Y)
            if (Xmax-Xmin)/dX+1 ~= numel(Y)
                error('X and Y arrays must be same length');
            end
            obj.Xmin = Xmin;
            obj.dX = dX;
            obj.Xmax = Xmax;
            obj.Y = Y;
        end
        
        % get all X values
        function x = get.X(obj)
            x = obj.Xmin:obj.dX:obj.Xmax;
        end
        
        % get minimum Y value
        function ym = get.Ymin(obj)
            ym = min(obj.Y(:));
        end
        
        % get maximum Y value
        function ym = get.Ymax(obj)
            ym = max(obj.Y(:));
        end
    end
end