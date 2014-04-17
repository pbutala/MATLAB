function oM = getCombinations(V,k)
if nargin == 1
    k = length(V);
end
if k>length(V)
    M = subFn2(V,k);
else
    M = repmat(V(:).',k,1);
    k2 = power(k,k);
    flags = ones(k2,1);
    if k > 1
        d = (k2-1)/(k-1);
        flags(1:d:end) = 0;
    end
    
    nck = nchoosek(V,k);
    [ITER,~] = size(nck);
    for i=1:ITER
        m = subFn(V(nck(i,:)));
        M(:,end+1:end+size(m,2)-k) = m(:,logical(flags));
    end
end
oM = removeDuplicates(M);
end

function M = subFn(V)
nV = length(V);
rnLn = power(nV,nV);
M = zeros(nV,power(nV,nV));
for iV = 1:nV
    stpLn = power(nV,nV-iV);
    iLn = 1;
    iter = 0;
    while(iLn <= rnLn)
        M(iV,iLn:iLn+stpLn-1) = V(rem(iter,nV)+1);
        iter = iter + 1;
        iLn = iLn+stpLn;
    end
end
end
function M = subFn2(V,k)
nV = length(V);
if k>nV
    rnLn = power(nV,k);
    M = zeros(k,rnLn);
    for iV = 1:k
        stpLn = power(nV,k-iV);
        iLn = 1;
        iter = 0;
        while(iLn <= rnLn)
            M(iV,iLn:iLn+stpLn-1) = V(rem(iter,nV)+1);
            iter = iter + 1;
            iLn = iLn+stpLn;
        end
    end
end
end

function M = removeDuplicates(V)
[~,nC] = size(V);
fOK = ones(1,nC);
for iC0 = 1:nC
    C0 = V(:,iC0);
    for iC1 = iC0+1:nC
        if fOK(iC1)
            C1 = V(:,iC1);
            fOK(iC1) = ~isequal(C0,C1);
        end
    end
end
M = V(:,logical(fOK));
end





































