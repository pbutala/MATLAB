classdef cSensor < handle
    properties
        location; % local
        dimension = cSize;
        responsivity = cCurve;
        filter = cCurve;
        type;
        color;
        % assuming each pixel is associated with a filter
        % if a global filter is applied, that can be part of cReceiver
    end
    
    properties(Dependent = true, SetAccess = private)
        % polygon ?
        % for now, treating pixels as rectangles
        pxR;
        pxT;
        pxL;
        pxB;
        edge;
    end
    
    methods
        % Constructor
        function obj = cSensor(Xs,Ys,Zs)
            if nargin~=0
                obj.location = cLocation(Xs,Ys,Zs);
            else
                obj.location = cLocation(0,0,0);
            end
            [m n] = size(Xs);
            obj.dimension.L(1:m,1:n) = 1e-3;
            obj.dimension.W(1:m,1:n) = 1e-3;
            obj.dimension.H(1:m,1:n) = 1e-3;
            obj.responsivity(1:m,1:n) = cCurve(200,1,1100,0.4*ones(1,901));
            obj.filter(1:m,1:n) = cCurve(200,1,1100,ones(1,901));
            obj.type = ones(m,n);
            obj.color = [1 1 1];
        end
        
        function FR = getResponsivityFiltered(obj)
            % gets filter.*responsivity for each pixel
            nP = obj.location.count;
            for iP = 1:1:nP
                fr = obj.filter(iP).*obj.responsivity(iP);
                FR(iP) = fr;
            end
        end
    end
    
    methods
        % set responsivity
        function set.responsivity(obj,value)
            if isa(value,'cCurve')
                obj.responsivity = value;
            end
        end
        
        % set filter
        function set.filter(obj,value)
            if isa(value,'cCurve')
                obj.filter = value;
            end
        end
        
        % set dimension
        function set.dimension(obj,value)
            if isa(value,'cSize')
                obj.dimension = value;
            end
        end
        
        function e = get.edge(obj)
            L = max(obj.location.X(:)) - min(obj.location.X(:))...
                + obj.dimension.L(1)/2 + obj.dimension.L(end)/2;
            W = max(obj.location.Y(:)) - min(obj.location.Y(:))...
                + obj.dimension.W(1)/2 + obj.dimension.W(end)/2;
            H = max(obj.dimension.H(:));
            e = cSize(L,W,H);
        end
        
        % get corner Rights
        function p = get.pxR(obj)
            p = obj.location.X + obj.dimension.L/2;
        end
        
        % get corner Top
        function p = get.pxT(obj)
            p = obj.location.Y + obj.dimension.W/2;
        end
        
        % get corner Left
        function p = get.pxL(obj)
            p = obj.location.X - obj.dimension.L/2;
        end
        
        % get corner Bottom
        function p = get.pxB(obj)
            p = obj.location.Y - obj.dimension.W/2;
        end
    end
end