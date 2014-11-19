function A = getSortedArray2(A)
% Y = getSortedArray(A)
% Uses merge-sort -> O(n)
% Returns A sorted in 'Ascending' (default mode) or 'Descending' order
A = mergesort(A,1,numel(A));
end

function X = mergesort(X,p,r)
if p<r
    q = floor((p+r)/2);
    X = mergesort(X,p,q);
    X = mergesort(X,q+1,r);
    X = merge(X,p,q,r);
end
end

function X = merge(X,p,q,r)
% merge iteration
Xl = X(p:q);
Nl = q-p+1;
Xr = X(q+1:r);
Nr = r-q;
i=1;j=1;
k=p;
while (i<=Nl) && (j<=Nr)
    if(Xl(i) < Xr(j))
        X(k) = Xl(i);
        i=i+1;
    else
        X(k) = Xr(j);
        j=j+1;
    end
    k=k+1;
end
if(i<=Nl)
    N = numel(Xl(i:end));
    X(k:k+N-1) = Xl(i:end);
else
    N = numel(Xr(j:end));
    X(k:k+N-1) = Xr(j:end);
end
end