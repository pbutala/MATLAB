% scrCSK
if exist('hWB','var') && ishandle(hWB)
    delete(hWB);
end
close all;
clearvars;
clc;

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP 4.Optimus
fSAVEALL = false;
fCLOSEALL = false;
fARCHIVE = false;
rng('default');

CHAROVERWRITE = '~';
STRPREFIX = '1_';
if(fARCHIVE)
    CHARIDXARCHIVE = '';           % ARCHIVE INDEX
else
    CHARIDXARCHIVE = CHAROVERWRITE; % OK TO OVERWRITE
end

%% CONSTANTS
LAMBDAMIN = 200; LAMBDADELTA = 1; LAMBDAMAX = 1100;
lambdas = LAMBDAMIN:LAMBDADELTA:LAMBDAMAX;

WBX = 50; WBY = 500; WBW = 275; WBH = 75;   % Wait Box X,,Y,WID,HGT
WBTITLE = 'Running CSK Simulation...'; % Wait Box title

%% M-CSK CONSTELLATION
% 4-CSK % xr = 0.734; yr = 0.265;  % xg = 0.402; yg = 0.597;
        % xb = 0.169; yb = 0.007;  % xw = 0.435; yw = 0.290;
%      00,    01,    10,    11     (g w b r)
x = [0.402; 0.435; 0.169; 0.734];
y = [0.597; 0.290; 0.007; 0.265];

% T = [xi xj xk; yi yj yk; 1 1 1];
% P = T\[xp;yp;1];

%% ranges

%% config
Yc = 400;        % Ilumination in CIE 1932 XYZ color model

% STATION
switch fSTATION
    % Results directory; Spource file; LED table dir;
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MatlabResults\14. CSK\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\scrCSK.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MatlabResults/14. CSK/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/scrCSK.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MatlabResults\\14. CSK\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSK.m';
    case 4
        ctDirRes = 'C:\\Users\\Pankil\\Documents\\MatlabResults\\14. CSK\\';
        ctFileCodeSrc = 'C:\\Users\\Pankil\\My Documents\\MATLAB\\Research\\Code V2\\Scripts\\scrCSK.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes STRPREFIX 'scrCSK' CHARIDXARCHIVE '.m']; % Script copy name
ctFileVars = [ctDirRes STRPREFIX 'datCSK' CHARIDXARCHIVE '.mat'];   % Data file name
if ~exist(ctDirRes,'dir')   % if data directory does NOT exist
    mkdir(ctDirRes);        % create data dir
end

%% logic
try
    %% prep start
    % Wait Bar to show progress
    hWB = waitbar(0,'Simulation: 0.00% done','Name',WBTITLE,...
        'Visible','Off',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    set(hWB,'Position',[WBX WBY WBW WBH],'Visible','On');
    setappdata(hWB,'canceling',0);
    LOOPCOUNT = 0; TOTALLOOPS = 2;
    TSTART = tic;
    
    %% compute
    CIDX = 3;
    xc = x(CIDX); yc = y(CIDX);
    Xc = Yc*xc/yc;
    Zc = Yc*(1-xc-yc)/yc;
    
    [Rc,Gc,Bc] = cCIE.xyz2rgb(Xc,Yc,Zc)
    sum = Rc + Gc + Bc
    [Xh,Yh,Zh] = cCIE.rgb2xyz(Rc,Gc,Bc)
    sumh = Xh + Yh + Zh
    
    %% prep stop
    % Update progress on wait bar
    LOOPCOUNT = LOOPCOUNT + 1;
    PROGRESS = LOOPCOUNT/TOTALLOOPS;
    TELAPSED = toc(TSTART);
    TREM = (TELAPSED/PROGRESS)-TELAPSED;
    waitbar(PROGRESS,hWB,sprintf('Simulation: %0.2f%% done...\nEstimated time remaining: %0.0f min',PROGRESS*100,TREM/60));
    if(getappdata(hWB,'canceling'))
        delete(hWB);
        error('Simulation aborted');
    end
    
    delete([ctDirRes '*' CHAROVERWRITE '.*']);
    if fSAVEALL                                                                 % SAVE script and data
        copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
        save(ctFileVars);                       % save workspace
    end
    delete(hWB);               
catch ex
    delete(hWB);
%     setpref('Internet','E_mail','pbutala@bu.edu');
%     setpref('Internet','SMTP_Server','smtp.bu.edu');
%     STREMAIL = ['Simulation ' STRPREFIX ' error!'];
%     sendmail('pankil.butala@gmail.com',STREMAIL);
    rethrow(ex);
end
% setpref('Internet','E_mail','pbutala@bu.edu');
% setpref('Internet','SMTP_Server','smtp.bu.edu');
% STREMAIL = ['Simulation ' STRPREFIX ' done. Starting plots.'];
% sendmail('pankil.butala@gmail.com',STREMAIL);
fprintf('--scrCSK Done--\n');































































