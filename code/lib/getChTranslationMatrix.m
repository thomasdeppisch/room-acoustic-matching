function T = getChTranslationMatrix(phiTrans, krTrans, chOrderIn, chOrderOut, indexing)
% get translation matrix in the CH domain using ACN channel ordering
% apply translation matrix on the left side, i.e.: pTrans = T * p;
%
% See e.g. Hahn, Spors, "Modal Bandwidth Reduction in Data-based Binaural
% Synthesis including Translatory Head-movements", DAGA.

% THIS USES/ASSUMES COMPLEX CHS!

if strcmpi(indexing, 'ACN')
    % ACN indexing
    mInTempVec = [-1:-1:-chOrderIn; 1:chOrderIn];
    mInIdxVec = [0; mInTempVec(:)];
    mOutTempVec = [-1:-1:-chOrderOut; 1:chOrderOut];
    mOutIdxVec = [0; mOutTempVec(:)];
else
    mInIdxVec = (-chOrderIn:chOrderIn)';
    mOutIdxVec = (-chOrderOut:chOrderOut)';
end

mInMinusMOutMtx = mInIdxVec.' - mOutIdxVec;
J = besselj(mInMinusMOutMtx, krTrans);
E = exp(1i * mInMinusMOutMtx * phiTrans);

T = (1i).^(mInMinusMOutMtx) .* J .* E;