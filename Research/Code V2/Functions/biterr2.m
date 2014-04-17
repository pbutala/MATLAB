%% function e = BITERR2(S1,S2)
% my custom function to incorporate biterr(...) functionality in R2011b
% version compared to R2013a

%% HEADER -----------------------------------------------------------------
% Author: Pankil Butala
% Email: pbutala@bu.edu 
% Institution: Multimedia Communications Laboratory,
%              Boston University.
% Generated: 16th August, 2013
% 
% Disclaimer: The file is provided 'as-is' to the user. It has not been
% tested for any bugs or inconsistent behavior. If you have any questions
% or bugs to report, please email the author at the specified email address
%
% Copyright (C) Pankil Butala 2013
% End Header --------------------------------------------------------------

function e = biterr2(S1,S2)
x = find(S1~=0 & S1~=1,1);
y = find(S2~=0 & S2~=1,1);
if ~isempty(x)|| ~isempty(y)
    error('input bit stream must contain only 0s or 1s');
end

e = sum(abs(S1(:) - S2(:)));
end