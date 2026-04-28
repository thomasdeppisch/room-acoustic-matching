function synthRir = synthIsoDiffRirChd(rirEstimate, dirSoundAziRad, fs, chOrder, targetRirLenSmp, freqRangeRt60Hz, preDirSoundLenMs, postDirSoundLenMs)
% synthesize a incoherent reverberation tail by using uncorrelated noise in
% the CH or SH domain
%
% This version applies same decay, same DRR, same energy to all channels
% (just based on omni channel of estimate)
arguments
    rirEstimate (:,:)
    dirSoundAziRad (:,:)
    fs (1,1)
    chOrder (1,1)
    targetRirLenSmp (1,1)
    freqRangeRt60Hz (1,2) = [100, 16000];
    preDirSoundLenMs (1,1) = 0.5;
    postDirSoundLenMs (1,1) = 1;
end

rirEstimate = rirEstimate(:,1); % only work with omni channel

numChs = 2*chOrder+1;
dirSound = getCH(chOrder,dirSoundAziRad,"real")';

% estimate RT60
octaveFactorRt60 = 1;
[rt60, fAxisIta, rirEstEndSmp] = getT20Ita(rirEstimate,fs,octaveFactorRt60,freqRangeRt60Hz);

% estimate DRR
maxTimeMsDirSound = 50;
numPeaksDirSoundToConsider = 1; % just take largest peak
[drrDb,dirSoundIdxSmp] = getDrr(rirEstimate(1:rirEstEndSmp),fs,preDirSoundLenMs,postDirSoundLenMs,freqRangeRt60Hz,maxTimeMsDirSound,numPeaksDirSoundToConsider);
drrLin = 10^(drrDb/10);

% figure
% plot(db(abs(rirEstimate)))
% hold on
% xline(dirSoundIdxSmp,'LineWidth',2)
% xline(rirEstEndSmp,'LineWidth',2)
% title('estimate')

% octave filter
filterOrder = 8;
freqRange = freqRangeRt60Hz;
fbBandwidth = '1 octave';
filtBank = octaveFilterBank(fbBandwidth, fs, 'FilterOrder', filterOrder, 'FrequencyRange', freqRange, 'OctaveRatioBase', 2);
filtBank2 = octaveFilterBank(fbBandwidth, fs, 'FilterOrder', filterOrder, 'FrequencyRange', freqRange, 'OctaveRatioBase', 2); % need this one with more channels

fcFiltBank = getCenterFrequencies(filtBank)';
numSubbands = length(fcFiltBank);
groupDelaySubbands = round(getGroupDelays(filtBank));

assert(all(fcFiltBank == fAxisIta), 'Filter bank frequencies must match the one from ITA!');

% split into direct and reverberant part
dirSoundStartIdx = max(1, dirSoundIdxSmp - round(preDirSoundLenMs/1000*fs) + 1);
dirSoundEndIdx = dirSoundIdxSmp + round(postDirSoundLenMs/1000*fs);

rirEstDirPart = rirEstimate(dirSoundStartIdx:dirSoundEndIdx);
rirEstRevPart = rirEstimate(dirSoundEndIdx+1:rirEstEndSmp);
numSamplesRevPart = targetRirLenSmp - 1; % we add the direct sound in the end

% zero pad, filter
numSamplesSubband = size(rirEstRevPart,1);
rirEstRevPartPadded = [rirEstRevPart;zeros(max(groupDelaySubbands),1)];
rirEstRevPartPaddedSubband = filtBank(rirEstRevPartPadded);

% create noise and decompose into subbands
noiseSig = randn(numSamplesRevPart, numChs);
numSamplesSubbandNoise = size(noiseSig,1);
noiseSigPadded = [noiseSig;zeros(max(groupDelaySubbands),numChs)];
noiseSigSubbandPadded = filtBank2(noiseSigPadded);

% reverberant energy estimation for each subband
noiseSubbandEnergyCorrected = zeros(numSamplesSubbandNoise, numSubbands, numChs);
tAxRevPart = (0:numSamplesRevPart-1).'/fs;

for ii = 1:numSubbands
    % pad to account for different group delays
    currSubbandRirRevAligned = squeeze(rirEstRevPartPaddedSubband(groupDelaySubbands(ii)+1:numSamplesSubband+groupDelaySubbands(ii),ii));
    revEnergyLin = sum(currSubbandRirRevAligned.^2);

    % same with noise
    noiseSubbandAligned = squeeze(noiseSigSubbandPadded(groupDelaySubbands(ii)+1:numSamplesSubbandNoise+groupDelaySubbands(ii),ii,:));

    % apply decay
    noiseSubbandDecaying = noiseSubbandAligned .* 10.^((-60./rt60(ii).*tAxRevPart) ./ 20);
    
    % adjust energy to energy of reverberant part of estimate
    noiseEnergyLin = mean(sum(noiseSubbandDecaying.^2));
    noiseSubbandEnergyCorrected(:,ii,:) = noiseSubbandDecaying .* sqrt(revEnergyLin ./ noiseEnergyLin); 
    % double check: mean(sum(squeeze(noiseSubbandEnergyCorrected(:,ii,:)).^2))

end

% back to time domain
shapedNoiseTimeDomain = squeeze(sum(noiseSubbandEnergyCorrected,2));
revEnergyResynth = mean(sum(shapedNoiseTimeDomain.^2));
dirEnergyTarget = drrLin * revEnergyResynth;

% add direct part
atfDirSoundOmniEnergy = dirSound(1).^2;
rirEstEarlyPart = dirSound * sqrt(dirEnergyTarget/atfDirSoundOmniEnergy);
% double check, this should be equal to drrDb: 10*log10(sum(rirEstEarlyPart(1).^2) / revEnergyResynth)

synthRir = [rirEstEarlyPart.'; shapedNoiseTimeDomain];
















