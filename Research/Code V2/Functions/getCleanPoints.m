%% function [Xo,Yo] = getCleanPoints(Xi,Yi,Zi,...,dmin)
% This function returns adjascent points on a plot that are greater than or
% equal to distance dmin. 

%% HEADER -----------------------------------------------------------------
% Author: Pankil Butala
% Email: pbutala@bu.edu 
% Institution: Multimedia Communications Laboratory,
%              Boston University.
% Generated: 25th October, 2013
% Modifications: 
% 
% Disclaimer: The file is provided 'as-is' to the user. It has not been
% tested for any bugs or inconsistent behavior. If you have any questions
% or bugs to report, please email the author at the specified email address
%
% Copyright (C) Pankil Butala 2013
% End Header --------------------------------------------------------------

function varargout = getCleanPoints(varargin)

NDIM = nargin-1;
LVECTI = length(varargin{1});
Pi = zeros(LVECTI,NDIM);
for i=1:NDIM
    V = varargin{i};
    Pi(:,i) = V(:);
end
d = varargin{NDIM+1};
Po = Pi(1,:);
Pt0 = Po;
for k = 2:LVECTI
   dP = sqrt(sum((Pi(k,:)-Pt0).^2,2));
   if dP >= d
       Po(end+1,:) = Pi(k,:);
       Pt0 = Pi(k,:);
   end
end
% if ~isequal(Pt0,Pi(end,:))
%     Po(end+1,:) = Pi(end,:);
%     Pt0 = Pi(k,:);
% end
for j=1:NDIM
    varargout{j} = Po(:,j);
end
end