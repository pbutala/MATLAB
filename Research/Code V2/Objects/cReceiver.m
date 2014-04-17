classdef cReceiver < handle
    properties
        location;
        orientation = cOrientation(0,0,0);
    end
    properties(Dependent = true, SetAccess = private)
        rxCount;
    end
    methods
        % Constructor
        function obj = cReceiver(Xs,Ys,Zs)
            obj.location = cLocation(Xs,Ys,Zs);
        end
        
        function n = get.rxCount(obj)
            n = obj.location.count;
        end
    end
    
    methods(Abstract)
        Isig = getSignal(obj,psd);
    end
end