classdef cOrientation
    properties(SetAccess = private)
        Z = 0;      % zenith angle
        A = 0;      % azimuth
        T = 0;      % tilt
%         N = eye(3); % basis of unit surface vector
    end
    
    methods
        % Constructor
        function obj = cOrientation(z,a,t)
            if nargin ~= 0
                if isequal(size(z),size(a),size(t))
                    obj.Z = z;
                    obj.A = a;
                    obj.T = t;
                else
                    error('Input arguments must have same dimensions')
                end
            end
        end
    end
end
