function filteredsignal = notchsignals(signal,fsamp)

% INPUTS:
%   signal - input signal (data structure supported by load_sig method);
%   fsamp - sampling frequency [Hz]
% OUTPUT
%   Y - (filtered) signal matrix, with each signal in separate row

frad = round(4/(fsamp/size(signal,2))); 
filteredsignal=zeros(size(signal));

for i=1:size(signal,1)                                    
    filteredsignal(i,:) = real(removelineinter(signal(i,:), frad, fsamp));  
end

function filteredsignal = removelineinter(signal,frad,fsamp)
% subtracts the content at the frequencies [Ind-r:Ind+r] from the signal y.
% OUTPUT
%   ynew - filtered signal

fsignal = fft(signal);
fcorrec = zeros(1,length(fsignal));  

window = fsamp;
tstamp = zeros(1,length(fsignal));
pos_tstamp = 1;

for i = 1:window:length(fsignal)-window
    
    mFreq = median(abs(fsignal(i+1:i+window)));
    sFreq = std(abs(fsignal(i+1:i+window)));
    tstamp2 = find(abs(fsignal(i+1:i+window)) > mFreq+5*sFreq);
    tstamp2 = tstamp2 + i;
    
    for j = -floor(frad/2):floor(frad/2)
        tstamp(pos_tstamp:pos_tstamp+length(tstamp2)-1) = tstamp2 + j;
        pos_tstamp = pos_tstamp + length(tstamp2);
    end

end

tstamp = round(tstamp(tstamp>0 & tstamp <= length(fsignal)/2+1));
for i = tstamp
    fcorrec(i) = fsignal(i);    
end

correc = length(fsignal) - floor(length(fsignal)/2) * 2;
fcorrec(length(fsignal):-1:ceil(length(fsignal)/2)+1) = conj(fcorrec(2:1:ceil(length(fsignal)/2)+1-correc));
filteredsignal = signal - ifft(fcorrec(1:length(fsignal)));

