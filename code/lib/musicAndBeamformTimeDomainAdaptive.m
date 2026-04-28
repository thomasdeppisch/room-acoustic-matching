function [sigBfd, medianDoaAziEleRad] = musicAndBeamformTimeDomainAdaptive(sig, stftBlockLenSmp, stftHopSizeSmp, ...
    doaObservationLenMs, arrayTransferFunctions, atfGridAziZenRad, params, trueDoaAziEleRad, plotFoundDoas)

if nargin < 7 || isempty(plotFoundDoas)
    plotFoundDoas = false;
end

stftWin = sqrt(hann(stftBlockLenSmp, 'periodic'));
fftLen = stftBlockLenSmp;
numPosFreqs = fftLen/2+1;
sigStft = stft(sig,stftWin,stftHopSizeSmp,fftLen,params.fs); % fftLen x numBlocks x numChannels
numBlocks = size(sigStft,2);
numMics = size(sigStft,3);

numFramesPerDoa = ceil((doaObservationLenMs/1000 * params.fs) / stftHopSizeSmp); % how many blocks should we observe for each doa estimate
if mod(numFramesPerDoa,2) == 0 % make odd number for symmetric centering
    numFramesPerDoa = numFramesPerDoa+1;
end
numAnalysisFramesOnEitherSide = (numFramesPerDoa-1)/2;

if ~isfield(params,'fftLen')
    params.fftLen = stftBlockLenSmp;
end

stftBfd = zeros(fftLen/2+1, numBlocks);
medianDoaAziEleRad = zeros(numBlocks,2);
for ii = 1:numBlocks
    % the frame we currently estimate the DOA for is ii, we center the
    % analysis around ii
    analysisStartFrame = max(ii-numAnalysisFramesOnEitherSide,1);
    analysisEndFrame = min(ii+numAnalysisFramesOnEitherSide,numBlocks);

    if ii < numAnalysisFramesOnEitherSide+1
        centerFrameIdx = ii;
    else
        centerFrameIdx = numAnalysisFramesOnEitherSide+1;
    end

    currSigStft = sigStft(1:numPosFreqs,analysisStartFrame:analysisEndFrame,:); % frames used in current estimation
    [currStftBfd, medianDoaAziEleRad(ii,:)] = musicAndBeamform(currSigStft, arrayTransferFunctions(1:numPosFreqs,:,:), atfGridAziZenRad, params, trueDoaAziEleRad, plotFoundDoas);
    stftBfd(:,ii) = currStftBfd(:,centerFrameIdx);
end

sigBfd = istft([stftBfd; conj(flipud(stftBfd(2:end-1,:)))], stftWin, stftHopSizeSmp, 'symmetric');


