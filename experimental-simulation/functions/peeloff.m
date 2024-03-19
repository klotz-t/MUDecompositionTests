function X = peeloff(X, spikes, fsamp, win)

windowl = round(win*fsamp);
waveform = zeros(1,windowl*2+1);
firings = zeros(1, size(X,2));
firings(spikes) = 1;
EMGtemp = zeros(size(X));

for l = 1:size(X,1)
    temp = cutMUAP(spikes,windowl,X(l,:));
    waveform = mean(temp,1);
    EMGtemp(l,:) = conv(firings,waveform,'same');
end 
    
X = X - EMGtemp;
