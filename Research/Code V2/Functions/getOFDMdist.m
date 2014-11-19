function [cdf, bins, lo, hi] = getOFDMdist(ofdmtyp,  nsc, syms, ofstDco, ofstAco, nbins, binfct, niter)

lo = realmax('double');
hi = -realmax('double');
if~exist('nbins','var')
    nbins = 1e4;
end
NBINS = nbins;
if(NBINS<1)
    error('Number of bins must be positive.');
end
if~exist('binfct','var')
    binfct = 1/1e3;
end
BINFCT = binfct; % 1%
if(BINFCT<=0)
    error('Bin factor must be positive.');
end
if~exist('niter','var')
    niter = power(2,10);
end
NITER = niter;
if(NITER<=0)
    error('Number of iterations must be positive.');
end
bins = 0:(1/(NBINS-1)):1;
hst = zeros(size(bins));


switch lower(ofdmtyp)
    case 'acoofdm'
        d = nsc/4;                % number of data carriers per ACOOFDM symbol
    case {'dcoofdm','dmt'}
        d = nsc/2 - 1;            % number of data carriers per DCOOFDM symbol
    otherwise
        error('OFDM type must be ''ACOOFDM'' or ''DCOOFDM'' or ''DMT''');
end
msc = numel(syms);
m = log2(msc);
iter = 1;
while iter < NITER
    txBits = randi([0 1],d,m);
    data = bin2decMat(txBits) + 1;
    % Generate OFDM symbol
    txSig = genOFDMsignal(... % Variable Arguments to the function
        'data',data,...
        'OFDMtype',ofdmtyp,...
        'N',nsc,...
        'Symbols',syms,...
        'ClipLow',0,...
        'ClipHigh',realmax('double'),...
        'OffsetDcoStddev', ofstDco,...
        'OffsetAcoStddev', ofstAco,...
        'ShowConst',false);
    
    sigL = min(txSig);
    sigH = max(txSig);
    if((sigH - hi) > (BINFCT*hi)) || ((sigL - lo) < (BINFCT*lo))
        clear pdf bins;
        tlo = min(lo, sigL);
        thi = max(hi, sigH);
        bins = tlo:((thi-tlo)/(NBINS-1)):thi;
        hst = zeros(size(bins));
        iter = 1;
    else
        hst = hst + hist(txSig, bins);
        iter = iter + 1;
    end
    lo = min(lo, sigL);
    hi = max(hi, sigH);
end
pdf = hst/sum(hst);
cdf = zeros(size(pdf));
cdf(1) = pdf(1);
for i=2:NBINS
    cdf(i) = cdf(i-1) + pdf(i);
end

if nargout == 0
    figure;
    subplot(2,1,1);
    plot(bins,pdf);
    title('PDF of OFDM signal values');
    xlabel('OFDM signal values');
    ylabel('pdf');
    grid on;
    
    subplot(2,1,2);
    plot(bins,cdf);
    title('CDF of OFDM signal values');
    xlabel('OFDM signal values');
    ylabel('cdf');
    grid on;
end
end