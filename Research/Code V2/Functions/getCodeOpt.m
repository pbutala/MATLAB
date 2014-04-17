function Copt = getCodeOpt(C)
% function Copt = getCodeOpt(C)
% returns distances between 
% -INPUT-
% -OUTPUT-

[~,nCODE] = size(C);
mD = zeros(nCODE,nCODE);        % init
for iC0 = 1:nCODE               % iter all codes
    C0 = C(:,iC0);              % select first
    for iC1 = iC0 : nCODE       % iter all 'next' codes
        C1 = C(:,iC1);          % select next
        V = C1-C0;              % calc vector from C0 to C1
        D = sqrt(sum(V.*V));    % calc distance from C0 to C1
        mD(iC0,iC1) = D;        % store distance in matrix
        mD(iC1,iC0) = D;        % store distance in matrix (across diag)
    end
end

% Copt = zeros(size(C));          % init optimal C
iCopt = zeros(1,nCODE);         % init index array for optimal C
[Rmax,Cmax] = find(mD == max(max(mD)),1); % find index of maximum distance codes     
iCopt(1) = Cmax;                % select Cth code
iCopt(2) = Rmax;                % select Rth code
for iC = 3:nCODE                            % iter to populate remaining indexes
    mseC = zeros(nCODE,1);                  % init MSE buffer
    for iC2 = 1:nCODE                       % iter thru all codes and chose one not selected in optimal yet
        if isempty(find(iCopt == iC2,1))    
            seC = 0;                        
            C2 = C(:,iC2);                  % C2 is new code
            for iC3 = 1:iC-1                % iter thru selected codes
                C3 = C(:,iCopt(iC3));       
                seC = seC + sum((C3-C2).^2);% calculate Sq Err
            end
            mseC(iC2) = seC/(iC-1);         % calculate MSE
        end
    end
    iCopt(iC) = find(mseC == max(mseC),1);  % select code with maximum MSE to optimal
end
Copt = C(:,iCopt);                          % convert optimal indexes to code and return
end