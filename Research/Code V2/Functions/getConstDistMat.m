function d = getConstDistMat(CNST)
    SYMNUM = size(CNST,2);
    d = zeros(SYMNUM);
    
    for i=1:SYMNUM
        c1 = CNST(:,i);
        for j=1:SYMNUM
            c2 = CNST(:,j);
            vc = c1-c2;
            d(i,j) = sqrt(sum(vc.*vc,1));
        end
    end
end