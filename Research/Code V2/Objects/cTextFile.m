classdef cTextFile < handle
    properties(SetAccess = immutable)
        path;
    end
    properties(SetAccess = private)
        handle;
    end
    
    methods
        function obj = cTextFile(path)
            obj.path = path;
        end
    end
end