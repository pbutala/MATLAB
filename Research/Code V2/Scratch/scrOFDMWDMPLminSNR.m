close all;
clc;
RNGSNROFST = RNGSNRDB - SNROFST;
for iTsd = 1:LENLEDWID                                                       % LOOP START LED SD
    for iT = 1:LENCCT
        for iFW = 1:LENFWHM                                                     % LOOP START FWHM
            for iNSC = 1:LENMODNSC                                                  % LOOP START NSC
                for iOf = 1:LENOFDMTYPES                                            % LOOP START OFDM types
                    for iTx = 1:NTX                                                 % LOOP START NTX
                        snrth = RNGSNROFST(find(BER(iTx,:,iTsd,iT,iFW,iNSC,iOf)<BERTH,1,'first'),iTsd,iT,iFW,iNSC,iOf);
                        if ~isempty(snrth)
                            SNRTH(iTx,iTsd,iT,iFW,iNSC,iOf) = snrth;
                        else
                            SNRTH(iTx,iTsd,iT,iFW,iNSC,iOf) = nan;
                        end
                    end % TX                                                        % LOOP STOP NTX
                    SNRTHBADTX(iTsd,iT,iFW,iNSC,iOf) = max(SNRTH(:,iTsd,iT,iFW,iNSC,iOf));
                end % OFDM                                                          % LOOP STOP OFDM types
            end % NSC                                                               % LOOP STOP NSC
        end % FWHM                                                               % LOOP STOP FWHM
    end % CCT                                                                   % LOOP STOP CCT
end % LED SD                                                                % LOOP STOP LED SD

% BEST
[V, I] = min(SNRTHBADTX(:));
[TSD,T,FW,~,OF] = ind2sub(size(SNRTHBADTX),I);
fprintf('Most Power Efficient Point\nSNR: %f dB\nCCT: %d K\nTXWID: %d nm\nFILTFWHM: %d nm\n%s\n\n',...
    V,RNGCCT(T),RNGLEDWID(TSD),RNGFWHM(FW),RNGOFDMTYPES{OF});

% WORST
[V, I] = max(SNRTHBADTX(:));
[TSD,T,FW,~,OF] = ind2sub(size(SNRTHBADTX),I);
fprintf('Least Power Efficient Point\nSNR: %f dB\nCCT: %d K\nTXWID: %d nm\nFILTFWHM: %d nm\n%s\n\n',...
    V,RNGCCT(T),RNGLEDWID(TSD),RNGFWHM(FW),RNGOFDMTYPES{OF});