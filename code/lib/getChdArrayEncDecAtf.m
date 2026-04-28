function [Dec,EncShort,EncIr,DecIrRaw] = getChdArrayEncDecAtf(atfs, atfGridHorAziRad, fs, chOrder, chDef, fftLen, encIrLen, regulConstEnc, doPlot)
% Calculate LS-optimal CHD encoder and decoder for arbitrary arrays.
% Based on Politis, Gamper, "﻿Comparing modeled and measurement-based 
% spherical harmonic encoding filters for spherical microphone arrays",
% WASPAA, 2017
%
% td, 2025

arguments
    atfs % numSamples x numMics x numAtfDirs
    atfGridHorAziRad
    fs
    chOrder
    chDef
    fftLen
    encIrLen
    regulConstEnc
    doPlot = false;
end

% Calculate array (decoder) matrix
numAtfDirs = size(atfs,3);
numMics = size(atfs,2);
numChs = 2*chOrder+1;
numPosFreqs = fftLen/2+1;

assert(numAtfDirs >= numChs)

fAx = linspace(0, fs/2, numPosFreqs)';
C = getCH(chOrder, atfGridHorAziRad, chDef);

Dec = zeros(numPosFreqs, numMics, numChs);
condD = zeros(numPosFreqs,1);
for ff = 1:numPosFreqs
    P = squeeze(atfs(ff,:,:)); % M x D
    Dec(ff,:,:) = P * C / (C' * C); % eq (23)
    condD(ff) = cond(squeeze(Dec(ff,:,:)));
end

% cond(C'*C)

if doPlot
    % figure
    % semilogx(fAx, condD)
    % grid on
    % title('condition number of array matrix D')
    
    figure
    semilogx(fAx,db(abs(Dec(:,:))))
    grid on
    title('decoding filts')
    xlim([20,20000])
end

%% use the estimated matrix as a basis for an encoder
I = eye(numChs, numChs);
Enc = zeros(numPosFreqs, numChs, numMics);
condE = zeros(numPosFreqs, 1);
for ff = 1:numPosFreqs
    currD = squeeze(Dec(ff,:,:));
    if numChs == 1
        currD = currD.';
    end
    Enc(ff,:,:) = I * currD' / (currD * currD' + regulConstEnc * eye(numMics)); % eq (19)
    condE(ff) = cond(squeeze(Enc(ff,:,:)));
end

if doPlot
    figure
    semilogx(fAx, condE)
    grid on
    title('condition number of regularized encoder E')
end
%% create time-domain encoding filters
EncIrRaw = zeros(fftLen, numChs, numMics);
DecIrRaw = zeros(fftLen, numMics, numChs);

for ii = 1:numMics
    if strcmp(chDef,'complex')
        EncIrRaw(:,:,ii) = ifft(getChFreqDomainConjugate(Enc(:,:,ii)));
        DecIrRaw(:,ii,:) = ifft(getChFreqDomainConjugate(squeeze(Dec(:,ii,:))));
    elseif strcmp(chDef,'real')
        EncIrRaw(:,:,ii) = ifft([Enc(:,:,ii); flipud(conj(Enc(2:end-1,:,ii)))], 'symmetric');
        DecIrRaw(:,ii,:) = ifft([Dec(:,ii,:); flipud(conj(Dec(2:end-1,ii,:)))], 'symmetric');
    end
end

EncIrRaw = circshift(EncIrRaw,fftLen/2,1);

EncIr = EncIrRaw(fftLen/2 - encIrLen/2 + 1 : fftLen/2 + encIrLen/2, :, :);
EncShort = fft(EncIr,fftLen);

% fadeInLenRefEncSmp = 200;
% fadeOutLenRefEncSmp = 200;
% fadeWinInRef = hann(2*fadeInLenRefEncSmp);
% fadeWinOutRef = hann(2*fadeOutLenRefEncSmp);
% EFilt(1:fadeInLenRefEncSmp,:,:) = EFilt(1:fadeInLenRefEncSmp,:,:) .* fadeWinInRef(1:fadeInLenRefEncSmp);
% EFilt(end-fadeOutLenRefEncSmp+1:end,:,:) = EFilt(end-fadeOutLenRefEncSmp+1:end,:,:) .* fadeWinOutRef(end-fadeOutLenRefEncSmp+1:end);

if doPlot
    figure
    subplot(211)
    plot(db(abs(EncIr(:,:))))
    grid on
    title('array encoding filters: time domain')
    
    subplot(212)
    semilogx(fAx,db(abs(EncShort(1:fftLen/2+1,:))))
    grid on
    xlim([20,20000])
    title('freq domain')
end