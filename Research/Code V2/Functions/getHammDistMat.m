function d = getHammDistMat(vctBitMap)
    SYMNUM = size(vctBitMap,1);
    d = zeros(SYMNUM);
    
    for i=1:SYMNUM
        for j=1:SYMNUM
            d(i,j) = biterr2(vctBitMap(i,:),vctBitMap(j,:));
        end
    end
end