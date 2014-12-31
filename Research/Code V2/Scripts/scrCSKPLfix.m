if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

CHAROVERWRITE = '';
STRPREFIX = 'CBCALL_';

%% set paths
ctFileCodeSrc = [mfilename('fullpath') '.m'];                           % get fullpath of current file
[ctScrDir,~,~] = fileparts(ctFileCodeSrc);                              % get scripts dir
cd(ctScrDir);                                                           % set scripts dir as pwd (reference)
addpath(genpath('..'));
ctDirRes = '..\..\..\..\MatlabResults\14. CSK\';
ctDirData = [ctDirRes STRPREFIX 'Data\'];
ctFileVars = [ctDirData STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];      % Data file name


scrCSKPL;