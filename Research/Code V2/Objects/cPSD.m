classdef cPSD < handle
    properties(SetAccess = immutable)
        npsd;
%         lmFlux;
%         rdFlux;
        % TODO: Add CIE.x and CIE.y values
    end
    
    properties(SetAccess = private)
        lmFlux;
        rdFlux;
    end
    
    methods
        % Constructor
        function obj = cPSD(xmin,delta,xmax,psd)
            if nargin ~= 0
                if min(psd)<0
                    error('Power Spectral Density cannot have negative values');
                end
                obj.rdFlux = sum(psd)*delta;
                if obj.rdFlux == 0
                    obj.npsd = cCurve(xmin,delta,xmax,psd);
                else
                    obj.npsd = cCurve(xmin,delta,xmax,psd/obj.rdFlux);
                end
                obj.npsd.Ylow = 0;
                obj.npsd.Yhigh = 0;
                V = getEyeSens(xmin,xmax,delta,1978);
                obj.lmFlux = 683*sum(psd.*V)*delta;
            end
        end
        
        function scaleOutputFlux(obj,scale)
            % scales the output flux of each channel by the scalar 'scale'
            % specified
            obj.lmFlux = scale*obj.lmFlux;
            obj.rdFlux = scale*obj.rdFlux;
        end
    end
    
    %% Operator Overloading
    methods
        % Overload '+' operator
        function res = plus(obj1,obj2)
            Xmin = min(obj1.npsd.Xmin, obj2.npsd.Xmin);
            dX = min(obj1.npsd.dX, obj2.npsd.dX);
            Xmax = max(obj1.npsd.Xmax, obj2.npsd.Xmax);
            Xs = Xmin:dX:Xmax;
            Y1s = obj1.npsd.getY(Xs);
            Y2s = obj2.npsd.getY(Xs);
            Ys = obj1.rdFlux*Y1s + obj2.rdFlux*Y2s;
            res = cPSD(Xmin,dX,Xmax,Ys);
        end
        
        % Overload '.*' operator
        function res = times(obj1,obj2) 
            if ~ismatrix(obj1) || ~ismatrix(obj2)
                error('Input arguments must be matrices');
            end
            
            if ~isequal(size(obj1),size(obj2))
                error('Input arguments must have same size');
            end
            
            if (isa(obj1,'cPSD') && isnumeric(obj2))
                obj = obj1;
                scl = obj2;
            elseif(isa(obj2,'cPSD') && isnumeric(obj1))
                obj = obj2;
                scl = obj1;
            else
                error('Input arguments must be of types numeric and cPSD');
            end
            [m,n] = size(obj);
            for i = 1:m
                for j = 1:n
                    mn = obj(i,j).npsd.Xmin;
                    d = obj(i,j).npsd.dX;
                    mx = obj(i,j).npsd.Xmax;
                    psd = (scl(i,j)*obj(i,j).rdFlux)*obj(i,j).npsd.Y;
                    res(i,j) = cPSD(mn,d,mx,psd);
                end
            end
        end
        
        % Overload './' operator
        function res = rdivide(obj1,obj2) 
            res = obj1.*(1/obj2);
        end
        
        % Overload '*' operator
        function res = mtimes(obj1,obj2) 
            if ~ismatrix(obj1) || ~ismatrix(obj2)
                error('Input arguments must be matrices');
            end
            
            if ~isequal(size(obj1,2),size(obj2,1))
                error('Matrix dimensions must agree [* n]*[n *]');
            end
                
            if (isa(obj1,'cPSD') && isnumeric(obj2)) ||...
               (isa(obj2,'cPSD') && isnumeric(obj1))
            else
                error('Input arguments must be of types numeric and cPSD');
            end
            
            r = size(obj1,1);
            n = size(obj1,2);
            c = size(obj2,2);
            for i=1:r
                row = obj1(i,:);
                for j=1:c
                    col = obj2(:,j);
                    vect = row.*col';
                    res(i,j) = vect(1);
                    for k=2:n
                        res(i,j) = res(i,j) + vect(k);
                    end
                end
            end
%             if(isa(obj1,'cPSD'))
%                 obj = obj1;
%                 scl = obj2;
%             else
%                 obj = obj2;
%                 scl = obj1;
%             end
%             
%             if(isscalar(scl))
%                 mn = obj.npsd.Xmin;
%                 d = obj.npsd.dX;
%                 mx = obj.npsd.Xmax;
%                 psd = (scl*obj.rdFlux)*obj.npsd.Y;
%                 res = cPSD(mn,d,mx,psd);
%             else
%                 error('mtimes not defined for objects of type %s',class(scl));
%             end
        end
    end
end
















