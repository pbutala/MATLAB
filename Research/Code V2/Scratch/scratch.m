close all;
clearvars;
clc;

res = 0:0.5:1;
resLEN = numel(res);
NCLR = 3;
IDX = ones(NCLR,1);

run = true;
while run
    disp(IDX');
    IDX(1) = IDX(1) + 1;
    for i=1:NCLR
        if IDX(i) > resLEN
            IDX(i) = 1;
            if i == NCLR
                run = false;
            else
                IDX(i+1) = IDX(i+1) + 1;
            end
        end
    end
end