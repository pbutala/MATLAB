classdef cLocation
    properties(SetAccess = private)
        X = 0;
        Y = 0;
        Z = 0;
    end
    
    properties(Dependent = true, SetAccess = private)
        count;
    end
    
    methods
        % Constructor
        function obj = cLocation(x,y,z)
            narginchk(0,3);
            if nargin ~= 0
                if isequal(size(x),size(y),size(z))
                    obj.X = x;
                    obj.Y = y;
                    obj.Z = z;
                else
                    error('Input arguments must have same dimensions')
                end
            end
        end
        
        % overload pperatoir '+'
        function L = plus(obj1,obj2)
            if isa(obj1,'cLocation') && isa(obj2,'cLocation')
                Xs = [obj1.X(:); obj2.X(:)];
                Ys = [obj1.Y(:); obj2.Y(:)];
                Zs = [obj1.Z(:); obj2.Z(:)];
                L = cLocation(Xs,Ys,Zs);
            elseif isa(obj1,'cLocation') && isempty(obj2)
                L = cLocation(obj1.X,obj1.Y,obj1.Z);
            elseif isempty(obj1) && isa(obj2,'cLocation')
                L = cLocation(obj2.X,obj2.Y,obj2.Z);
            else
                error('Both input arguments to plus must be ''cLocation''');
            end
        end
        % return number of locations
        function n = get.count(obj)
            n = numel(obj.X);
        end
    end
end