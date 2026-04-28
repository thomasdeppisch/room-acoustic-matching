function bnCirc = circModalCoeffs(N, kr, arrayType, harmonicsDef)
% this is an equivalent to sphModalCoeffs for circular/equatorial arrays

bnSph = sphModalCoeffs(N, kr, arrayType);
Nnm = getNnm(N,pi/2,harmonicsDef);

bnCirc = zeros(size(bnSph, 1), 2*N+1);
for mm = -N : N
    if mm==0
        mIdxAcn = 1;
    else
        mIdxAcn = 2*abs(mm) + double(mm>0);
    end
    %fprintf('m=%i, idx=%i \n', mm, mIdxAcn)
    
    for nn = abs(mm) : N
        bnCirc(:, mIdxAcn) = bnCirc(:, mIdxAcn) + bnSph(:, nn+1) .* Nnm(nn^2+nn+mm+1).^2;
    end
end