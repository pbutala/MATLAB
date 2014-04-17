% scratch CIE plot
close all;
clearvars;
clc;

figure;
imshow('figCIEXYZ.png');
axis tight;
axis on;
% x = 0.195 + 0.5*(0.835-0.195)/0.8;
% y = 0.165 + 0.5*(0.915-0.165)/0.9;
% annotation('ellipse',[x y 0.01 0.01],'EdgeColor','k','LineWidth',2);
Xlim = get(gca,'Xlim');
Ylim = get(gca,'Ylim');

ah = annotation('ellipse','EdgeColor','k','LineWidth',2);
set(ah,'parent',gca);
set(ah,'position',[32 12 10 10]);

% L = 380:5:780;
% X = zeros(size(L));
% Y = zeros(size(L));
% figure;
% set(gca,'units','normalized');
% hold on;
% obs = cCIE;
% for iL = 1:numel(L)
%     psd = cCurve(L(iL),1,L(iL),1);
%     [X(iL) Y(iL)] = obs.getCoordinates(psd);
% end
% plot([X X(1)],[Y Y(1)]);
% axis([0 1 0 1]);
% axis equal;

% % some data with an awkward axis
% plot(30:40,rand(1,11));
% % create a text annotation
% ah=annotation('textbox');
% % ...and position it on the current axis
% set(ah,'parent',gca);
% set(ah,'position',[31 .1 3 .2]);
