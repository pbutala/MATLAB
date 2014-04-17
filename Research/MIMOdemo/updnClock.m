function Y = updnClock(X,xFs,cFs,show)
[USF,DSF] = rat(cFs/xFs);

% upsample
uX = upsample(X,USF);
% take fourier transform of upsampled signal
uL = length(uX);
uN = 2^nextpow2(uL);
uFT = USF*fft(uX,uN)/uL;
uf = (xFs*USF/2)*linspace(0,1,uN/2+1);
uFTM = abs(uFT);

% create filter
flt = zeros(uN,1);
flt(uf<=xFs/2) = 1;
flt(uf>xFs*USF-xFs/2) = 1;

% filter upsampled signal
iFT = uFT.*flt;
iFTM = abs(iFT);
iL = uL;
iN = 2^nextpow2(iL);
iX = iL*ifft(iFT,iN,1,'symmetric');
ufX = iX(1:uL);

% Downsample
Y = ufX(1:DSF:end);

if nargin ~= 4
    show = false;
end
if show == true
    figure('Name',sprintf('Change Sample Rate, up:%d dn:%d',USF,DSF),'NumberTitle','OFF');
    % signal
    subplot(3,4,1);
    plot(X);
    axis tight;
    title(sprintf('X, sample rate = %0.2f Msps',xFs*1e-6'));
    ax(1) = subplot(3,4,5);
    ax(2) = subplot(3,4,9);
    plotFreq(X,xFs,'single',ax);
    
    % upsampled
    subplot(3,4,2);
    plot(uX);
    axis tight;
    title(sprintf('X^{up}, USF = %d',USF));
    ax(1) = subplot(3,4,6);
    ax(2) = subplot(3,4,10);
    plotFreq(uX,xFs*USF,'single',ax);
    
    % upsampled filtered
    subplot(3,4,3);
    plot(ufX);
    axis tight;
    title(sprintf('X^{ip}, sample rate = %0.2f Msps',xFs*USF*1e-6));
    ax(1) = subplot(3,4,7);
    ax(2) = subplot(3,4,11);
    plotFreq(ufX,xFs*USF,'single',ax);
    
    % Y
    subplot(3,4,4);
    plot(Y);
    axis tight;
    title(sprintf('Y, sample rate = %0.2f Msps',cFs*1e-6'));
    ax(1) = subplot(3,4,8);
    ax(2) = subplot(3,4,12);
    plotFreq(Y,cFs,'single',ax);
    
%     % take fourier transform of signal
%     sL = length(X);
%     sN = 2^nextpow2(sL);
%     sFT = fft(X,sN)/sL;
%     sf = xFs/2*linspace(0,1,sN/2+1);
%     sFTM = abs(sFT);
%     
%     % take fourier transform of signal
%     yL = length(Y);
%     yN = 2^nextpow2(yL);
%     yFT = fft(Y,yN)/yL;
%     yf = cFs/2*linspace(0,1,yN/2+1);
%     yFTM = abs(yFT);
%     
%     figure('Name',sprintf('Change Sample Rate, up:%d dn:%d',USF,DSF),'NumberTitle','OFF');
%     subplot(2,4,1);
%     plot(X);
%     axis tight;
%     title(sprintf('X, sample rate = %0.2f Msps',xFs*1e-6'));
%     
%     subplot(2,4,5);
%     plot(sf*1e-6,2*sFTM(1:sN/2+1));
%     axis tight;
%     title('Single sided amplitude spectrum of X');
%     xlabel('Frequency (MHz)');
%     ylabel('|X_f|');
%     
%     subplot(2,4,2);
%     plot(uX);
%     axis tight;
%     title(sprintf('X^{up}, USF = %d',USF));
%     
%     subplot(2,4,6);
%     ax = plotyy(uf*1e-6,2*uFTM(1:uN/2+1),uf*1e-6,flt(1:uN/2+1));
%     axis(ax, 'tight');
%     title('Single sided amplitude spectrum of upsampled X');
%     xlabel('Frequency (MHz)');
%     ylabel('|X_f^{up}|');
%     
%     subplot(2,4,3);
%     plot(ufX);
%     axis tight;
%     title(sprintf('X^{ip}, sample rate = %0.2f Msps',xFs*USF*1e-6));
%     
%     subplot(2,4,7);
%     plot(uf*1e-6,2*iFTM(1:iN/2+1));
%     title('Single sided amplitude spectrum of filtered upsampled signal');
%     xlabel('Frequency (MHz)');
%     ylabel('|X_f^{ip}|');
%     axis tight;
%     
%     subplot(2,4,4);
%     plot(Y);
%     axis tight;
%     title(sprintf('Y, sample rate = %0.2f Msps',cFs*1e-6'));
%     
%     subplot(2,4,8);
%     plot(yf*1e-6,2*yFTM(1:yN/2+1));
%     axis tight;
%     title('Single sided amplitude spectrum of Y');
%     xlabel('Frequency (MHz)');
%     ylabel('|Y_f|');
end