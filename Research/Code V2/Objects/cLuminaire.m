classdef cLuminaire < handle
    properties 
        orientation = cOrientation(pi,0,0);
        dimension = cSize(20e-2,20e-2,5e-2);
        order = 1;
        % transmission;
        % driver ? 
        % bandwidth ?
    end
    
    properties (SetAccess = immutable)
        location;
        channels; % color channels
    end
    
    properties (Dependent = true, SetAccess = private)
        lmFlux;
        rdFlux;
        lmFluxClr;
        rdFluxClr;
        color;
        chCount;
        lmCount;
    end
    methods
        % Constructor
        function obj = cLuminaire(clrChnls,Xs,Ys,Zs)
            if isa(clrChnls,'cPSD')
                obj.channels = clrChnls;
            else
                error('Argument ''clrChnls'' must be of type cPSD');
            end
            if nargin == 2
                Ys = Xs.Y;
                Zs = Xs.Z;
                Xs = Xs.X;
            end
            [m,n] = size(Xs);
            obj.location = cLocation(Xs,Ys,Zs);
            obj.orientation = cOrientation(pi,0,0);
            D = 20e-2*ones(m,n);
            H = 5e-2*ones(m,n);
            obj.dimension = cSize(D,D,H);
        end
        % Cannot calculate the free space gain here because luminaire does
        % not have information about room reflections
        
        function scaleOutputFlux(obj,scale)
            % scales the output flux of each channel by the scalar 'scale'
            % specified
            nC = obj.chCount;
            for iC = 1:nC
                obj.channels(iC).scaleOutputFlux(scale);
            end
        end
    end
    
    %% Get/Set accessors
    methods
        % get luminous flux for each color in matrix
        function f = get.lmFluxClr(obj)
            nC = numel(obj.channels);
            f = zeros(size(obj.channels));
            for iC = 1:nC
                f(iC) = obj.channels(iC).lmFlux;
            end
        end
        
        % get radiant flux for each color in matrix
        function f = get.rdFluxClr(obj)
            nC = numel(obj.channels);
            f = zeros(size(obj.channels));
            for iC = 1:nC
                f(iC) = obj.channels(iC).rdFlux;
            end
        end
        
        % get total luminous flux
        function f = get.lmFlux(obj)
            fC = obj.lmFluxClr;
            f = sum(fC(:));
        end
        
        % get total radiant flux
        function f = get.rdFlux(obj)
            fC = obj.rdFluxClr;
            f = sum(fC(:));
        end
        
        % get color output of the luminaire
        function clr = get.color(obj)
            nC = numel(obj.channels);
            clr = cPSD(0,1,0,0);
            for iC = 1:nC
                clr = clr + obj.channels(iC);
            end
        end
        
        % get number of channels for each luminaire
        function n = get.chCount(obj)
            n = numel(obj.channels);
        end
        
        % get number of luminaires of this type
        function n = get.lmCount(obj)
            n = obj.location.count;
        end
    end
end
    