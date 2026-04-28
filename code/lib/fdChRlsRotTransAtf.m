function [ir, h, rmsDiffDbSmooth, isConverged, convTimeSec, allIrs] = fdChRlsRotTransAtf(mcSigIn, refSigIn, fs, chOrder, chDef, ...
    arrayMatrix, rotationAngleAziRad, transAziRad, transRadius, regulConst, ...
    forgettingFact, blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly, convThreshDb)
% Frequency-domain RLS to CH domain for moving arrays, uses array transfer
% functions instead of analytic model.

arguments
    mcSigIn (:,:)
    refSigIn (:,1)
    fs (1,1)
    chOrder (1,1)
    chDef {mustBeMember(chDef,["real","complex"])}
    arrayMatrix (:,:,:) % numFreqs x numMics x numChs
    rotationAngleAziRad (:,1)
    transAziRad (:,1)
    transRadius (:,1)
    regulConst  
    forgettingFact = [] % [] .. cumulative average, otherwise exponential smoothing where 1 means no forgetting, 0 means all forgetting
    blockLenSmp (1,1) = 0.5*fs
    hopSizeSmp (1,1) = 2048
    filterLenSmp (1,1) = blockLenSmp/2
    win (:,1) = sqrt(hann(blockLenSmp, 'periodic')); % hann works better than rect
    stopEarly = false;
    convThreshDb = -20;
end

if size(mcSigIn,1) ~= size(refSigIn, 1)
    warning('Shortening signal length.')
    minLen = min(size(mcSigIn,1), size(refSigIn, 1));
    mcSigIn = mcSigIn(1:minLen, :);
    refSigIn = refSigIn(1:minLen);
end

assert(length(mcSigIn) == length(refSigIn) && ...
    length(mcSigIn) == length(rotationAngleAziRad) && ...
    length(mcSigIn) == length(transAziRad) && ...
    length(mcSigIn) == length(transRadius))

fftLen = blockLenSmp;
if mod(fftLen, 2) == 1 % enforce even FFT length
    fftLen = fftLen + 1;
end

numFreqs = fftLen;
numPosFreqs = fftLen/2+1;
overlapSmp = blockLenSmp-hopSizeSmp;
numBlocks = floor((size(mcSigIn,1) - overlapSmp) / hopSizeSmp);
numMics = size(mcSigIn, 2);

if isscalar(regulConst)
    regulConst = repmat(regulConst, numPosFreqs, 1); % use frequency dependent regularization
end

chOrderCombi = chOrder; % TODO: we might want this as input parameter
numChannelsChs = 2*chOrderCombi+1;

c = 343;
fAx = linspace(0,fs/2,numPosFreqs)';
k = 2*pi*fAx/c;

saveAllIrs = false;
if nargout > 3
    saveAllIrs = true;
    allIrs = zeros(2*(numPosFreqs-1),numChannelsChs,numBlocks);
end

h = zeros(numPosFreqs, numChannelsChs);
Rxx = zeros(numPosFreqs, numChannelsChs, numChannelsChs);
rxd = zeros(numPosFreqs, numChannelsChs);
%pDiff = zeros(numPosFreqs, numMics); % diffuse field response
%pDiffTOnly = zeros(numPosFreqs, numChannelsChs); % diffuse field response modification due to translation

if stopEarly
    rmsDiffDbSmooth = 0;
    hOld = 0;
    expAvgFact = 0.99;
end

isConverged = false;

warning('off','MATLAB:nearlySingularMatrix'); % turn off warning avoid spamming the terminal
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:illConditionedMatrix');

%figure

for bb = 1:numBlocks
    blockIdx = 1 + ((bb-1) * hopSizeSmp) + (0:blockLenSmp-1);
    mcSigBlock = fft(mcSigIn(blockIdx, :) .* win, fftLen);
    refSigBlock = fft(refSigIn(blockIdx) .* win, fftLen);

    % rotation matrix
    circMeanRad = circularMeanAziRad(rotationAngleAziRad(blockIdx));
    rotMtx = getChRotMtx(chOrder, circMeanRad, chDef);

    % average translation parameters within block
    transAziMeanRad = circularMeanAziRad(transAziRad(blockIdx));
    transRadiusMean = mean(transRadius(blockIdx));

    for ff = 1:numPosFreqs
        currK = 2*pi*fAx(ff)/c;
        currD = squeeze(arrayMatrix(ff,:,:)); % numMics x numChs

        if numChannelsChs == 1
            currD = currD.';
        end
        pCurr = mcSigBlock(ff,:).';

        T = getChTranslationMatrix(transAziMeanRad, currK*transRadiusMean, chOrderCombi, chOrder, 'ACN'); 

        X = currD * rotMtx * T * refSigBlock(ff);
        % Rxx(ff,:,:) = forgettingFact .* squeeze(Rxx(ff,:,:)) + X' * X;
        % rxd(ff,:) = forgettingFact .* rxd(ff,:).' + X' * pCurr; 

        %pDiffCurr = real(diag(currD * rotMtx * (T * T') * rotMtx' * currD'));

        if isempty(forgettingFact) || forgettingFact == 1
            % use a cumulative average!
            Rxx(ff,:,:) = (bb-1)/bb .* squeeze(Rxx(ff,:,:)) + 1/bb * (X' * X);
            rxd(ff,:) = (bb-1)/bb .* rxd(ff,:).' + 1/bb * X' * pCurr; 
        else
            % use exponential smoothing
            Rxx(ff,:,:) = forgettingFact .* squeeze(Rxx(ff,:,:)) + (1-forgettingFact) * (X' * X);
            rxd(ff,:) = forgettingFact .* rxd(ff,:).' + (1-forgettingFact) * X' * pCurr;

            % or classical RLS:
            % Rxx(ff,:,:) = forgettingFact .* squeeze(Rxx(ff,:,:)) + (X' * X);
            % rxd(ff,:) = forgettingFact .* rxd(ff,:).' + X' * pCurr;
            
        end

        % if the regulConst has 2 values, we treat the first one as out-of-bound
        % and the second one as in-bound value
        h(ff,:) = (squeeze(Rxx(ff,:,:)) + eye(numChannelsChs)*regulConst(ff)) \ rxd(ff,:).';
    end

    if stopEarly
        rmsDiffDbSmooth = expAvgFact * rmsDiffDbSmooth + (1-expAvgFact) * db(rms(abs(h - hOld), 'all'));   

        % check parseval for early stopping
        % plot(bb,rmsDiffDbSmooth,'*')
        % hold on

        hOld = h;
        if rmsDiffDbSmooth < convThreshDb
            disp(['RLS converged after ' num2str(bb) '/' num2str(numBlocks) ' blocks'])
            isConverged = true;
            break
        end
    end

    if saveAllIrs
        hDouble = getChFreqDomainConjugate(h);
        allIrs(:,:,bb) = ifft(hDouble);
    end
end

warning('on','MATLAB:nearlySingularMatrix'); % turn warning on again
warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:illConditionedMatrix');

hDouble = getChFreqDomainConjugate(h);
ir = ifft(hDouble);
ir = ir(1:filterLenSmp, :);

if saveAllIrs
    allIrs = allIrs(1:filterLenSmp, :, :);
end

convTimeSec = (blockLenSmp + (bb-1) * hopSizeSmp) / fs;


% figure
% semilogx(fAx,db(abs(h.')))
% grid on
% xlim([20,20000])
