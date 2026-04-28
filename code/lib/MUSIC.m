function [doa, steeringVecDoa] = MUSIC(signalFrame, numDOAsToEstimate, arrayTransferFunctions, gridXYZ)
% signalFrame .. numSamples x numChannels
% arrayTransferFunctions .. numDirs x numChannels
% gridXYZ .. numDirs x 3

    [Uc,~,~] = svd(signalFrame.');
    Un = Uc(:, numDOAsToEstimate+1:end); % noise subspace
    
    VnA = Un'*arrayTransferFunctions.';
    P_music = 1./sum(conj(VnA).*VnA).';

    peaksIdx = peakFind2d(P_music, [gridXYZ(:,1), gridXYZ(:,2), gridXYZ(:,3)], numDOAsToEstimate);
    doa = [gridXYZ(peaksIdx,1), gridXYZ(peaksIdx,2), gridXYZ(peaksIdx,3)].';
    steeringVecDoa = arrayTransferFunctions(peaksIdx,:).';
end