close all;
clearvars;
clc;
% scrFigNameChange
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\11. SISOFDM all\';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/11. SISOFDM all/';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\11. SISOFDM all\\';
    otherwise
        error('Station not defined');
end
files = dir(ctDirRes);
NF = numel(files);
for i=1:NF
    fSRC = [ctDirRes files(i).name];
    I = find(fSRC=='~');
    if ~isempty(I) 
        fDEST = [fSRC(1:I-1),fSRC(I+1:end)];
        copyfile(fSRC,fDEST);
        delete(fSRC);
    end
end
% delete([ctDirRes '*' '~' '.*']);