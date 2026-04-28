function [dereverbSig, dereverbStft, origStft] = adaptiveGwpeDerev(sig, params)
% Generalized Weighted Prediction Error (GWPE)
% This function only allows the identity matrix approach for covariance
% estimation, but it is about 4 times faster than the "full" implementation.
% This is a block-based implementation, able to adapt to changes in the
% transfer functions.

blockLenStft = params.blockLenStft;
hopSizeStft = params.hopSizeStft;
numBlocksCovEst = params.numBlocksCovEst; % within these macro blocks, the covariance will be estimated
fftLen = params.fftLen;
predOrders = params.predOrders;
predDelay = params.predDelay;
numIterations = params.numIterations;
covSmoothingLen = params.covSmoothingLen;
regulWeight = params.regulWeight;
adaptiveCovInit = params.adaptiveCovInit;

% calculate which samples to save in overlap save
overlapSmpCov = round(params.overlapFactCov * numBlocksCovEst);
if mod(overlapSmpCov,2) == 1
    overlapSmpCov = overlapSmpCov + 1; % ensure even overlap
end
numDiscardSmp = overlapSmpCov/2;
numSaveSmp = numBlocksCovEst - 2*numDiscardSmp;
hopSizeCov = numBlocksCovEst - overlapSmpCov;

if size(predOrders,1) == 1
    predOrders = repmat(predOrders,[fftLen/2+1,1]);
end

win = sqrt(hann(blockLenStft, 'periodic'));
%ffts = stft(sig,'window',win,'overlapLength',blockLen-hopsize,'FFTLength',fftLen,'FrequencyRange','onesided');
ffts = stft(sig,win,hopSizeStft,fftLen,params.fs);
ffts = ffts(1:fftLen/2+1,:,:);
origStft = ffts;

ffts = permute(ffts,[3 2 1]);

[numChannels, numBlocks, numFreqs] = size(ffts);
x = ffts;

parfor ff = 1:numFreqs
%for ff = 1:numFreqs
    predOrder = predOrders(ff);
    y = squeeze(ffts(:, :, ff));

    numSubBlocks = numBlocks - predOrder + 1;
    numMacroBlocks = floor((numBlocks - predOrder - predDelay + 1 - overlapSmpCov) / hopSizeCov);

    yBuff = zeros(numChannels, predOrder, numSubBlocks);
    for cc = 1:numChannels
        yBuff(cc, :, :)= flipud(buffer(y(cc,:), predOrder, predOrder-1, 'nodelay')); % flip time, y(t), ..., y(t-K+1)
    end
    yBuff = reshape(yBuff, numChannels*predOrder, numSubBlocks); % mono-frequent signal over time in blocks with hopsize 1, flipped in time (new to old samples)

    % these temp variables are unfortunately needed for pafor
    tempX = x(:,:,ff); % this is the working copy
    tempXClean = x(:,:,ff); % this is a clean copy
    tempXOut = zeros(size(tempX)); % this will be written to the output
    tempGy = zeros(numChannels, numBlocks);
    g = zeros(predOrder*numChannels,numChannels);
    for bb = 1:numMacroBlocks
        delStartIdx = (bb-1) * hopSizeCov + 1; % indices of the "delayed" blocks
        delEndIdx = min(delStartIdx + numBlocksCovEst - 1, numSubBlocks);
        nonDelStartIdx = delStartIdx + predDelay + predOrder - 1;
        nonDelEndIdx = min(nonDelStartIdx+numBlocksCovEst-1, numBlocks);
        if bb == 1 % save full first half
            saveSmpIdxNonDel = nonDelStartIdx : nonDelEndIdx-numDiscardSmp;
        elseif bb == numMacroBlocks % save full second half
            saveSmpIdxNonDel = nonDelStartIdx+numDiscardSmp : nonDelEndIdx;
        else % only save central part (overlap save)
            saveSmpIdxNonDel = nonDelStartIdx+numDiscardSmp : nonDelStartIdx+numDiscardSmp+numSaveSmp-1;
        end
        %saveSmpIdxDel = saveSmpIdxNonDel - predOrder - predDelay + 1;

        if adaptiveCovInit
            % This uses the old signal variance estimate to initialize the
            % new estimation. It can lead to very aggressive filtering
            % resulting in musical noise.
            tempX(:, nonDelStartIdx:nonDelEndIdx) = y(:, nonDelStartIdx:nonDelEndIdx) - g.' * yBuff(:, delStartIdx:delEndIdx); % first update x with the old filter coefficients
        else
            % This resets the signal variance in each block, giving a
            % result closer to the original GWPE
            tempX(:, nonDelStartIdx:nonDelEndIdx) = tempXClean(:, nonDelStartIdx:nonDelEndIdx);
        end

        for iter=1:numIterations
            lambda = mean(abs(tempX(:, nonDelStartIdx:nonDelEndIdx)).^2, 1); % covariance estimate via identity matrix
            if covSmoothingLen > 1
                lambda = movmean(lambda,covSmoothingLen,2);
            end
    
            %invlambda = 1 ./ max(lambda, regulWeight); 
            invLambda = 1 ./ (lambda + regulWeight); % this regularization makes it equivalent to the "full" implementation
            
            yDel = yBuff(:, delStartIdx:delEndIdx);
            tmp = (conj(yDel) .* invLambda);
            R = 1/numBlocksCovEst * tmp * yDel.';
            r = 1/numBlocksCovEst * tmp * y(:, nonDelStartIdx:nonDelEndIdx).';
            g = (R + eye(numChannels*predOrder)*regulWeight*trace(R)/(numChannels*predOrder)) \ r;
            tempGy(:, delStartIdx:delEndIdx) = g.' * yDel;
            tempX(:, nonDelStartIdx:nonDelEndIdx) = y(:, nonDelStartIdx:nonDelEndIdx) - tempGy(:, delStartIdx:delEndIdx); % x is iteratively optimized, bufferedSubchannels are not!
        end
        tempXOut(:, saveSmpIdxNonDel) = tempX(:, saveSmpIdxNonDel);
    end

    x(:,:,ff) = tempXOut;
end

dereverbStft = permute(x,[3 2 1]);
%dereverbSig = istft(dereverbStft,'window',win,'overlapLength',blockLen-hopsize,'FFTLength',fftLen,'FrequencyRange','onesided');
%reverbSig = istft(permute(Gy,[3 2 1]),'window',win,'overlapLength',blockLen-hopsize,'FFTLength',fftLen,'FrequencyRange','onesided');

dereverbSig = istft([dereverbStft; conj(flipud(dereverbStft(2:end-1,:,:)))], win, hopSizeStft, 'symmetric');

end

