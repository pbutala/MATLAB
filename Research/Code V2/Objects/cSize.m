classdef cSize
    properties
        L = 5e-3;
        W = 5e-3;
        H = 2e-3;
    end
    
    properties (Dependent = true, SetAccess = private)
        A
        V
    end
    
    methods 
        % Constructor
        function obj = cSize(l,w,h)
            if nargin ~= 0
                if isequal(size(l),size(w),size(h))
                    obj.L = l;
                    obj.W = w;
                    obj.H = h;
                else
                    error('Input arguments must have same dimensions')
                end
            end
        end
        
        % returns area of object L*W
        function A = get.A(obj)
            A = obj.L.*obj.W;
        end
        % returns volume of object L*W*H
        function A = get.V(obj)
            A = obj.L.*obj.W.*obj.H;
        end
        
        % overload pperatoir '+'
        function D = plus(obj1,obj2)
            if isa(obj1,'cSize') && isa(obj2,'cSize')
                Ls = [obj1.L(:); obj2.L(:)];
                Ws = [obj1.W(:); obj2.W(:)];
                Hs = [obj1.H(:); obj2.H(:)];
                D = cSize(Ls,Ws,Hs);
            elseif isa(obj1,'cSize') && isempty(obj2)
                D = cSize(obj1.L,obj1.W,obj1.H);
            elseif isempty(obj1) && isa(obj2,'cSize')
                D = cSize(obj2.L,obj2.W,obj2.H);
            else
                error('Both input arguments to plus must be ''cSize''');
            end
        end
    end
end