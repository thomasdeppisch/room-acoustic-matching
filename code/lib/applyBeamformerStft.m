function [stftBfd, bfWeights] = applyBeamformerStft(sigStft, steeringVecs, bfType, regulConst, noisePsdMtx)
arguments
    sigStft (:,:,:) % numFrequencies (singleSided) x numBlocks x numChannels
    steeringVecs (:,:) % numChannels x numFrequencies or numChannels x 1
    bfType {mustBeMember(bfType,{'delayAndSum','delayAndSumNormd','MPDR','MVDR','wMPDR'})} = 'delayAndSumNormd';
    regulConst = 1e-6;
    noisePsdMtx = []; % numChannels x numChannels x numFrequencies
end

[numFrequencies, numBlocks, numChannels] = size(sigStft);
if size(steeringVecs,2) == 1
    steeringVecs = repmat(steeringVecs, [1 numFrequencies]);
end

stftBfd = zeros(numFrequencies, numBlocks); % median doa
bfWeights = zeros(numChannels, numFrequencies);

if nargin > 4 && ~isempty(noisePsdMtx)
    if strcmp(noisePsdMtx, 'identity')
        noisePsdMtx = repmat(eye(numChannels), [1 1 numFrequencies]);
    end
end

%regulConstMpdr = regulConst * max(sum(abs(sigStft).^2,[2,3]));
for ii = 1:numFrequencies
    if strcmp(bfType, 'delayAndSum') % Van Trees eq. 2.32
        bfWeights(:,ii) = steeringVecs(:,ii);
    elseif strcmp(bfType, 'delayAndSumNormd') % aka matched filter
        bfWeights(:,ii) = steeringVecs(:,ii) / (steeringVecs(:,ii)' * steeringVecs(:,ii) + regulConst * max(sum(abs(steeringVecs).^2,2)));
    elseif strcmp(bfType, 'MPDR') % Van Trees
        sig = squeeze(sigStft(ii,:,:));
        numBlocksCovEst = size(sig,1);

        Rxx = 1/numBlocksCovEst * (sig' * sig) + 1e-10 * eye(numChannels);
        bfWeights(:,ii) = (Rxx \ steeringVecs(:,ii)) / ...
            (steeringVecs(:,ii)' / Rxx * steeringVecs(:,ii) + regulConst);
    elseif strcmp(bfType, 'MVDR') % Van Trees, ch. 6
        bfWeights(:,ii) = (noisePsdMtx(:,:,ii) \ steeringVecs(:,ii)) / ...
            (steeringVecs(:,ii)' / noisePsdMtx(:,:,ii) * steeringVecs(:,ii) + regulConst);
    elseif strcmp(bfType, 'wMPDR') % See e.g. Boedekker et al., "Jointly optimal dereverberation and beamforming"
        sig = squeeze(sigStft(ii,:,:)).';
        lambda = mean(abs(sig).^2, 1);
        numBlocksCovEst = size(sig,1);

        for iter=1:2
            RMpdr = 1/numBlocksCovEst * (sig ./ (lambda+1e-4) * sig');

            bfWeights(:,ii) = (RMpdr \ steeringVecs(:,ii)) / ...
                (steeringVecs(:,ii)' / RMpdr * steeringVecs(:,ii) + regulConst);

            lambda = abs(bfWeights(:,ii)' * sig).^2;
        end
    end

    stftBfd(ii,:) = bfWeights(:,ii)' * squeeze(sigStft(ii,:,:)).';
end

