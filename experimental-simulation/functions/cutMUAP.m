%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracts cosecutive MUAPs out of signal Y and stores them row-wise in the
% output variable MUAPs
% From Farina's lab

% INPUTS:
%   - MUPulses      triggering positions (in samples) for rectangular window used in extraction of MUAPs (MU firring patterns);
%   - len           radius of rectangular window (window length = 2*len+1)
%   - Y             single signal channel (raw vector containing a single channel of a recorded signals)
%
% OUTPUTS:
%   - MUAPs         row-wise matrix of extracted MUAPs (aligned signal intervals of length 2*len+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MUAPs=cutMUAP(MUPulses,len,Y)

    while ~isempty(MUPulses) && MUPulses(end)+2*len>length(Y)
        MUPulses = MUPulses(1:end-1);
    end
    c=length(MUPulses);
    edgeLen = round(len/2);
    tmp = gausswin(2*edgeLen)';
    win = [tmp(1:edgeLen) ones(1,2*len-2*edgeLen+1) tmp(edgeLen+1:end)];
    MUAPs=zeros(c,1+2*len);
    for k=1:c
        MUAPs(k,1:1+2*len)= win .*[ zeros(1,max(MUPulses(k)-len,1)-(MUPulses(k)-len)) Y(max(MUPulses(k)-len,1):min(MUPulses(k)+len,length(Y))) zeros(1,MUPulses(k)+len-min(MUPulses(k)+len,length(Y)))];
    end
end