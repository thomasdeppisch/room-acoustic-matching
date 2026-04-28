function [stftBfd, medianDoaAziEleRad, bfWeights, atfsMedian] = musicAndBeamform(sigStft, arrayTransferFunctions, atfGridAziZenRad, params, trueDoaAziEleRad, plotFoundDoas)
arguments
    sigStft % numPosFreqs x numBlocks x numMics
    arrayTransferFunctions % numPosFreqs x numMics x numAtfDirs
    atfGridAziZenRad % numAtfDirs x 2
    params
    trueDoaAziEleRad
    plotFoundDoas
end

if nargin < 6 || isempty(plotFoundDoas)
    plotFoundDoas = false;
end

arrayTransferFunctions = permute(arrayTransferFunctions, [2,3,1]);
[numChannels, numPwDirs, numPosFrequencies] = size(arrayTransferFunctions);

assert(numPosFrequencies == size(sigStft,1)) % atfs and STFT need same frequency axis

numDoasToEstimate = 1;
[atfGridXYZ(:,1), atfGridXYZ(:,2), atfGridXYZ(:,3)] = sph2cart(atfGridAziZenRad(:,1), pi/2-atfGridAziZenRad(:,2), ones(numPwDirs,1));

if nargin < 5 || isempty(trueDoaAziEleRad)
    doas = zeros(numPosFrequencies,3);
    fAx = linspace(0,params.fs/2,params.fftLen/2+1).';
    freqRangeLogicals = fAx>params.musicFreqRange(1) & fAx<params.musicFreqRange(2);
    freqBinIdx = (1:numPosFrequencies)';
    freqRangeIdx = freqBinIdx(freqRangeLogicals);
    for ii = freqRangeIdx'
        currAtfs = arrayTransferFunctions(:,:,ii);
        doas(ii,:) = MUSIC(squeeze(sigStft(ii,:,:)), numDoasToEstimate, currAtfs.', atfGridXYZ);
    end
    doas = doas(freqRangeIdx,:);

    medianDoaCart = median(doas);
    medianDoaCart = medianDoaCart./sqrt(sum(medianDoaCart.^2,2));
    [medianDoaAziEleRad(1), medianDoaAziEleRad(2)] = cart2sph(medianDoaCart(1), medianDoaCart(2), medianDoaCart(3));

    % doaAziRad = cart2sph(doas(:,1), doas(:,2), doas(:,3)) * 180/pi;
    
    if plotFoundDoas
        figure
        [xSph,ySph,zSph] = sphere;
        surf(xSph,ySph,zSph,'FaceColor','none')
        hold on
        plot3(doas(:,1),doas(:,2),doas(:,3),'*')
        plot3(medianDoaCart(1),medianDoaCart(2),medianDoaCart(3),'r*')
        grid on
        axis equal
        xlabel('x')
        ylabel('y')
        title('MUSIC: found DOAs')
    end

else
    warning('musicAndBeamform: true DOA provided, skipping DOA estimation')
    medianDoaAziEleRad = trueDoaAziEleRad;
end


[medianDoaCart(1), medianDoaCart(2), medianDoaCart(3)] = sph2cart(medianDoaAziEleRad(1), medianDoaAziEleRad(2), 1);
[~,atfIdxClosestToMedianDoa] = min(sqrt(sum((medianDoaCart - atfGridXYZ).^2, 2)));
% atfGridAziZenRad(atfIdxClosestToMedianDoa,1) * 180/pi

atfsMedian = squeeze(arrayTransferFunctions(:,atfIdxClosestToMedianDoa,:));

if strcmp(params.bfType,'MVDR') && isfield(params, 'noisePsdMtx')
    [stftBfd, bfWeights] = applyBeamformerStft(sigStft(1:numPosFrequencies,:,:), atfsMedian, params.bfType, params.regulConst, params.noisePsdMtx);
else
    [stftBfd, bfWeights] = applyBeamformerStft(sigStft(1:numPosFrequencies,:,:), atfsMedian, params.bfType, params.regulConst);
end

