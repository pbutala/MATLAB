classdef cImagingOptics < handle
    properties
        fov = pi/3;
        f = 1e-2;
        fN = 1;
        T = 1;
    end
    properties (Dependent = true, SetAccess = private)
        D
        Area;
    end
    methods
        % Constructor
        function obj = cImagingOptics(f,fN)
            if nargin ~= 0
                obj.f = f;
                obj.fN = fN;
            end
        end
        % get aperture diameter
        function diam = get.D(obj)
            diam = obj.f/obj.fN;
        end
        
        % get aperture area (pi*D*D/4)
        function area = get.Area(obj)
            area = pi*(obj.D^2)/4.0;
        end
    end
end