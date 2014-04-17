% Hany's file
clc
close all
clear all

% FLAGS
fSTATION = 1;   % 1.PHO445 2.ENGGRID 3.LAPTOP
fSAVEALL = true;
fCLOSEFIGS = false;

% STATION
switch fSTATION
    case 1
        ctDirRes = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Results\5. OSM\';
        ctFileCodeSrc = '\\ad\eng\users\p\b\pbutala\My Documents\MATLAB\Research\Code V2\Scripts\OSM.m';
    case 2
        ctDirRes = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Results/5. OSM/';
        ctFileCodeSrc = '/home/pbutala/My Documents/MATLAB/Research/Code V2/Scripts/OSM.m';
    case 3
        ctDirRes = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Results\\5. OSM\\';
        ctFileCodeSrc = 'C:\\Users\\pbutala\\Documents\\MATLAB\\Research\\Code V2\\Scripts\\OSM.m';
    otherwise
        error('Station not defined');
end
ctFileCodeDest = [ctDirRes 'scrOSM.m'];
ctFileVars = [ctDirRes 'OSM.mat'];      % file to store workspace

m=2; % bits per symbols
M = 2^m; % # of symbols
alphabit = modem.pammod('M', 2^m); %PAM modulation

Nt=4; % # 0f Txs 
Nr = 4; % # of Rxs
snr =0:5:30; % signal-to-noise ratio
loop =1;

for jj =snr
    err = 0; % Detected Errors
    iter = 1000; % # of Iterations
    for ii = 1:iter
        % bit_stream = randint(1, log2(M*Nt)); % generate bit stream e.g. 4 bits per symbol for 4 Txs
        bit_stream = randi([0,1], [1,log2(M*Nt)]);
        sym_index = bin2dec(int2str(bit_stream(1:log2(M)))); % Symbol index
        led_index = bin2dec(int2str(bit_stream(log2(M)+1:log2(Nt*M))))+1; % LED index
        tx = zeros(Nt,1);
        tx(led_index)=modulate(alphabit,sym_index); % assign Symbol to LED 
        
        %%%%%%%%Imaging Rx Channel Martix%%%%%%%%%%%%%%%%
        Htx = abs(real(tx))';
        Htxindex=find(tx~=0);
        Htx(Htxindex)=1;
        H= zeros(Nt,Nr);       
        H(Htxindex,:)=Htx; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%Add AWGN Channel%%%%%%%%%%%%%%%%%%%%%%%%%%
        y= awgn((H*tx),jj,'measured');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%Optimal Spatial Modulation Detector%%%%%%%%%%%%%%
        for i=1:Nt
            for j=1:M
                g_jq = H(:,i)*modulate(alphabit,j-1);
                min_2(i,j) = norm(g_jq, 'fro')^2 - 2*real(y'*g_jq); %This vector contains the errors to all antennas.
            end
        end
        [tmp,idx] =  min( min_2(:)) ;
        [idx1,idx2] = find(min_2==tmp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%Estimated bits from LED index%%%%%%%%%%%%%%%%%%%
        tmp1 = dec2bin(idx1-1,log2(Nt));
        
        for uu = 1:log2(Nt)
            led_bits(uu)= bin2dec(tmp1(uu));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%Estimated bits from Symbol%%%%%%%%%%%%%%%%%%%%%%
        tmp2 = dec2bin(idx2-1,log2(M));
        
        for xx = 1:log2(M)
            sym_bits(xx)= bin2dec(tmp2(xx));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rx_bits = [sym_bits led_bits];

        err = err+biterr(rx_bits,bit_stream);
    end
    
        bit_err(loop) = err/((log2(M*Nt))*iter)
        loop = loop+1
end

    % delete old files
if fSAVEALL
    delete([ctDirRes '*.png']);
    delete([ctDirRes '*.mat']);
    delete([ctDirRes '*.fig']);
end

% save data
if fSAVEALL
    copyfile(ctFileCodeSrc,ctFileCodeDest); % save script
    save(ctFileVars);                       % save workspace
end

figure
hold on
semilogy(snr,bit_err,'-bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5)
xlabel('SNR (dB)')
ylabel('Bit Error Ratio')
% title('Imaging Rx OSM 4x4')
grid
hold off

if fSAVEALL
    fname = [ctDirRes ' BER vs SNR.fig'];
    saveas(gcf,fname,'fig');
    fname = [ctDirRes ' BER vs SNR.png'];
    saveas(gcf,fname,'png');
end
if fCLOSEFIGS
    close;
end