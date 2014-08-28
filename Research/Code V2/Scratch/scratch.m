close all;
clearvars;
clc;

A = 1:100;
IDStart = 1;
while numel(A) > 1
    IDSrvr = IDStart:2:numel(A);
    if (IDSrvr(end) == numel(A))
        IDStart = 2;
    else
        IDStart = 1;
    end
    A = A(IDSrvr)
end