function t0Smp = findDirSoundPeakIdx(ir,fs,maxTimeMs,numPeaks)
% find dir sound as earliest of largest N peaks
arguments
    ir
    fs (1,1)
    maxTimeMs (1,1) = 20
    numPeaks (1,1) = 5
end

maxTimeSmp = round(maxTimeMs/1000*fs);

[peakVals, peakIdx] = findpeaks(abs(ir(1:maxTimeSmp)), 'MinPeakDistance', 10, 'SortStr', 'descend', 'NPeaks', numPeaks);
[~,firstPeakSubIdx] = min(peakIdx);
t0Smp = peakIdx(firstPeakSubIdx);

