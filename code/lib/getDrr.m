function [drrDb,t0Smp] = getDrr(ir,fs,preDirSoundLenMs,postDirSoundLenMs,freqRangeHz,maxTimeMsDirSound,numPeaksDirSoundToConsider)
arguments
    ir
    fs
    preDirSoundLenMs
    postDirSoundLenMs
    freqRangeHz
    maxTimeMsDirSound = 20
    numPeaksDirSoundToConsider = 3;
end
numChannels = size(ir,2);

if nargin < 4 || isempty(freqRangeHz)
    freqRangeHz = [50, 16000];
end

preDirSoundLenSmp = round(preDirSoundLenMs/1000*fs);
postDirSoundLenSmp = round(postDirSoundLenMs/1000*fs);

[bBp,aBp] = butter(4, freqRangeHz/fs*2);
irFilt = filtfilt(bBp,aBp,ir);

drrDb = zeros(numChannels,1);
t0Smp = zeros(numChannels,1);

for chIdx = 1:numChannels
    t0Smp(chIdx) = findDirSoundPeakIdx(irFilt(:,chIdx), fs, maxTimeMsDirSound, numPeaksDirSoundToConsider);
    [~,rirEndSmp] = getRtItaBroadband(irFilt(:,chIdx), fs, freqRangeHz);

    dirPart = irFilt(max(1, t0Smp(chIdx)-preDirSoundLenSmp):t0Smp(chIdx)+postDirSoundLenSmp, chIdx);
    revPart = irFilt(t0Smp(chIdx)+postDirSoundLenSmp+1:rirEndSmp, chIdx);

    drrDb(chIdx) = 10*log10(sum(abs(dirPart).^2) ./ sum(abs(revPart).^2));

end