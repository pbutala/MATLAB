classdef cMM < cCBC
    % class to handle Metameric Modulation with CBCs
    properties(SetAccess = immutable)
        M;
        N;
    end % properties - immutable
    
    properties(Constant, GetAccess = private)
%         CBT = [6,2,0;6,1,0;5,2,0;5,1,0;4,2,0;4,1,0;3,2,0;3,1,0;2,1,0];
        CBT = [6,2,0;6,1,0;5,2,0;5,1,0;4,2,0;4,1,0;3,2,0;3,1,0]; % exclude CBC9.
    end % properties - immutable, private
    
    methods
        function obj = cMM(M,N)
            % obj = cMM(M) is class constructor
            if ~exist('M','var')
                M = 4;
            end
            if ~exist('N','var')
                N = 5;
            end
            switch M
                case{2;4;8} 
                    obj.M = M;
                otherwise
                    error('%d-MM is not supported',M);
            end
            obj.N = N;
        end
        
        function CBCs = getCBCs(obj)
            switch obj.M
                case 2
                    CBCs = obj.getCBCsM2();
                case 4
                    CBCs = obj.getCBCsM4();
                case 8
                    CBCs = obj.getCBCsM8();
            end
        end % getColorBands
    end % methods
    
    methods(Access=private)
        % M=2 gets all combinations of CBCs for different N
        function CBCs = getCBCsM2(obj)
            CBCs = [];
            for iC1 = 1:size(obj.CBT,1);
                for iC2 = iC1+1:size(obj.CBT,1)
                    if(numel(unique(obj.CBT([iC1,iC2],:))) == obj.N)
                        CBCs(end+1,:) = [iC1,iC2];
                    end
                end % iC2
            end % iC1
            if isempty(CBCs)
                error('%d-MM with N=%d is not supported',obj.M,obj.N);
            end
        end % getCBCsM2
        
        % M=4 gets all combinations of CBCs for different N
        function CBCs = getCBCsM4(obj)
            CBCs = [];
            for iC1 = 1:size(obj.CBT,1);
                for iC2 = iC1+1:size(obj.CBT,1)
                    for iC3 = iC2+1:size(obj.CBT,1);
                        for iC4 = iC3+1:size(obj.CBT,1)
                            if(numel(unique(obj.CBT([iC1,iC2,iC3,iC4],:))) == obj.N)
                                CBCs(end+1,:) = [iC1,iC2,iC3,iC4];
                            end
                        end % iC4
                    end % iC3
                end % iC3
            end % iC1
            if isempty(CBCs)
                error('%d-MM with N=%d is not supported',obj.M,obj.N);
            end
        end % getCBCsM4
        
        function CBCs = getCBCsM8(obj)
            CBCs = [];
            for iC1 = 1:size(obj.CBT,1);
                for iC2 = iC1+1:size(obj.CBT,1)
                    for iC3 = iC2+1:size(obj.CBT,1);
                        for iC4 = iC3+1:size(obj.CBT,1)
                            for iC5 = iC4+1:size(obj.CBT,1);
                                for iC6 = iC5+1:size(obj.CBT,1)
                                    for iC7 = iC6+1:size(obj.CBT,1);
                                        for iC8 = iC7+1:size(obj.CBT,1)
                                            if(numel(unique(obj.CBT([iC1,iC2,iC3,iC4,iC5,iC6,iC7,iC8],:))) == obj.N)
                                                CBCs(end+1,:) = [iC1,iC2,iC3,iC4,iC5,iC6,iC7,iC8];
                                            end
                                        end % iC8
                                    end % iC7
                                end % iC6
                            end % iC5
                        end % iC4
                    end % iC3
                end % iC2
            end % iC1
            if isempty(CBCs)
                error('%d-MM with N=%d is not supported',obj.M,obj.N);
            end
        end % getCBCsM8
    end % methods - private
end













































