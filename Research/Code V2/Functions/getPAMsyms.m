%% Syms = getPAMsyms(M,D,P) returns PAM symbols 
% getPAMsyms(M,D,P) plots PAM symbols on a scatterplot.
%
% This function returns normalized symbols for M-PAM where M is an integral
% exponent of 2. D (default = 1) specifies the delta between PAM symbols. P
% (Default = 0) is the polarity flag. -1:negative|0:bipolar|1:positive.
% If no output variable specified, the function displays
% the constellation on a scatterplot.

%% HEADER -----------------------------------------------------------------
% Author: Pankil Butala
% Email: pbutala@bu.edu
% Institution: Multimedia Communications Laboratory,
%              Boston University.
% Last Modified: 12th August, 2013
%
% Disclaimer: The file is provided 'as-is' to the user. It has not been
% tested for any bugs or inconsistent behavior. If you have any questions
% or bugs to report, please email the author at the specified email address
%
% Copyright (C) Pankil Butala 2013
% End Header --------------------------------------------------------------

function Syms = getPAMsyms(M,D,P)
bps = log2(M);      % bits per symbol

% Sanity check for integral exponent (bps must be integer)
if (rem(bps,1) ~= 0) || (M==0)  
%     error('''M'' must be an integral exponent of 2');
end
if ~exist('D','var')
    D = 1;                      % delta between PAM symbols
end
if ~exist('P','var')
    P = 0;                      % polarity of PAM symbols -1:negative, 0:bipolar, 1:positive
end

switch P
    case -1
        Sym0 = -(M-1)*D;        % smallest symbol
        Sym1 =  0;              % largest symbol
    case 0
        Sym0 = -(M-1)*D/2;       % smallest symbol
        Sym1 = (M-1)*D/2;        % largest symbol
    case 1
        Sym0 = 0;              % smallest symbol
        Sym1 = (M-1)*D;        % largest symbol
    otherwise
        error('Polarity P must be either -1,0 or 1');
end
Syms = (Sym0:D:Sym1)';          % Generate symbols

if nargout == 0                 % if no output argument specified, display constellation
    figure;                     % Generate axes for scatterplot
    scatter(Syms,zeros(length(Syms),1));             % Display the constellation on current axes
    set(gca,'XTick',Syms);       % Set X tick values
    grid on;                    % Show the grid
    axis([Sym0-D Sym1+D -1 1]); % Scale axes for pleasant display
    tStr = sprintf('%d-PAM constellation diagram',M);   % Generate title
    title(tStr);                % Show title
end