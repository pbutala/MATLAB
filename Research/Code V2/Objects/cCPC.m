classdef cCPC < handle
    properties
        fov = pi/2;
        index = 1.5;
    end
    properties (Dependent = true, SetAccess = private)
        gain;
    end
    
    methods
        % Constructor
        function obj = cCPC(index,fov)
            if nargin ~= 0
                obj.index = index;
                obj.fov = fov;
            end
        end
        % get CPC gain
        function G = get.gain(obj)
            G = (obj.index*obj.index)/(sin(obj.fov)*sin(obj.fov));
        end
    end
end