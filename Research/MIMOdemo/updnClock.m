function Y = updnClock(X,xFs,cFs,fltTyp,show)
if ~exist('fltTyp','var')
    fltTyp = 'IdealRect';
end
if ~exist('show','var')
    show = false;
end
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
switch lower(fltTyp)
    case 'idealrect'
        flt(uf<=xFs/2) = 1;
        flt(uf>xFs*USF-xFs/2) = 1;
    case 'raisedcosine'
        flt(uf<=xFs) = sinc(uf(uf<=xFs)/xFs);
    otherwise
        error('FLTTYP must be either ''IDEALRECT'' or ''RAISEDCOSINE''.');
end

% filter upsampled signal
iFT = uFT.*flt;
iFTM = abs(iFT);
iL = uL;
iN = 2^nextpow2(iL);
iX = iL*ifft(iFT,iN,1,'symmetric');
ufX = iX(1:uL);

% Downsample
Y = ufX(1:DSF:end);

% if nargin ~= 4
%     show = false;
% end
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
end