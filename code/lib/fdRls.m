function [ir, isConverged, convTimeSec, allIrs] = fdRls(inMc, inBfd, fs, forgettingFact, blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly)
% Frequency-domain RS for blind system identification after Meyer-Kahlen, 
% Schlecht, "Blind Directional Room Impulse Response Parameterization from 
% Relative Transfer Functions", 2023.
% Thomas Deppisch, 2023

arguments
    inMc (:,:)
    inBfd (:,1)
    fs (1,1)
    forgettingFact (1,1) = 0.99;
    blockLenSmp (1,1) = 0.5*fs;
    hopSizeSmp (1,1) = 2048;
    filterLenSmp (1,1) = blockLenSmp/2;
    win (:,1) = ones(blockLenSmp, 1); %win = sqrt(hann(blockLenSmp, 'periodic'));
    stopEarly = false;
end

if size(inMc,1) ~= size(inBfd, 1)
    warning('Shortening signal length.')
    minLen = min(size(inMc,1), size(inBfd, 1));
    inMc = inMc(1:minLen, :);
    inBfd = inBfd(1:minLen, :);
end

fftLen = blockLenSmp;
if mod(fftLen, 2) == 1 % enforce even FFT length
    fftLen = fftLen + 1;
end

numFreqs = fftLen/2+1;
overlapSmp = blockLenSmp-hopSizeSmp;
numBlocks = floor((size(inMc,1) - overlapSmp) / hopSizeSmp);
numMics = size(inMc, 2);

if stopEarly
    normDiffDbSmooth = 0;
    expAvgFact = 0.99;
    convThreshDb = -20;
    oldHNorm = 0;
end

isConverged = false;

saveAllIrs = false;
if nargout > 3
    saveAllIrs = true;
    allIrs = zeros(2*(numFreqs-1),numMics,numBlocks);
end

h = zeros(numFreqs, numMics);
phi = zeros(numFreqs, 1);
for bb = 1:numBlocks
    blockIdx = 1 + ((bb-1) * hopSizeSmp) + (0:blockLenSmp-1);
    inMcBlock = fft(inMc(blockIdx, :) .* win, fftLen);
    inBfdBlock = fft(inBfd(blockIdx, :) .* win, fftLen);
    inMcBlock = inMcBlock(1:numFreqs, :);
    inBfdBlock = inBfdBlock(1:numFreqs, :);

    epsilon = inMcBlock - h .* inBfdBlock;
    phi = forgettingFact * phi + conj(inBfdBlock) .* inBfdBlock;
    h = h + 1./phi .* conj(inBfdBlock) .* epsilon;

    if stopEarly
        % check parseval for early stopping
        currHNorm = norm(h,'fro');
        normDiffDbSmooth = expAvgFact * normDiffDbSmooth + (1-expAvgFact) * db(abs(currHNorm - oldHNorm));   
        oldHNorm = currHNorm;
        if normDiffDbSmooth < convThreshDb
            disp(['RLS converged after ' num2str(bb) '/' num2str(numBlocks) ' blocks'])
            isConverged = true;
            break
        end
    end

    if saveAllIrs
        allIrs(:,:,bb) = ifft([h; flipud(conj(h(2: end-1,: )))]);
    end
end

ir = ifft([h; flipud(conj(h(2: end-1,: )))]);

% shiftLenSmp = 100;
% ir = circshift(ir, shiftLenSmp);
% 
% fadeInLen = 60;
% winIn = hann(2*fadeInLen);
% winIn = winIn(1:fadeInLen);
% 
% ir(1:fadeInLen, :) = ir(1:fadeInLen,:) .* winIn;
ir = ir(1:filterLenSmp, :);

convTimeSec = (blockLenSmp + (bb-1) * hopSizeSmp) / fs;

