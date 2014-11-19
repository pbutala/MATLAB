% plot freq spectrum of the signal
function plotFreq(X,clk,side,ax,style)
if ~exist('side','var')
    side = 'single';
end
if ~exist('style','var')
    style = 'b-';
end
% take fourier transform of signal
L = length(X);
N = 2^nextpow2(L);
FT = fft(X,N)/L;
FTM = abs(FT);
FTA = angle(FT);

fr = clk/2*linspace(0,1,floor(N/2)+1);
switch(lower(side))
    case 'single'
        f = fr;
        mag = FTM(1:floor(N/2)+1);
        ang = FTA(1:floor(N/2)+1);
    case 'whole'
        f = [-fr(end-1:-1:2) fr];
        mag = circshift(FTM(:),floor(N/2)-1);
        ang = circshift(FTA(:),floor(N/2)-1);
    otherwise
        error('Side must be either ''Single'' or ''Whole''');
end

if ~exist('ax','var')
    figure;
    ax(1) = subplot(2,1,1);
    ax(2) = subplot(2,1,2);
end
axes(ax(1));
plot(f*1e-6,mag,style);
xlabel('Frequency (MHz)');
ylabel('|X_f|');
title('Magnitude');
axis tight;

if length(ax) > 1
    axes(ax(2));
    plot(f*1e-6,180*ang/pi,style);
    xlabel('Frequency (MHz)');
    ylabel('\theta_f^0');
    title('Phase (degrees)');
    axis tight;
end