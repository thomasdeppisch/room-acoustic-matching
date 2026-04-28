% This script demonstrates the RIR estimation and resynthesis procedure
% from Deppisch et al., "Identification and Matching of Room Acoustics With
% Moving Head-Worn Microphone Arrays", 2025.

clear all
close all

addpath(genpath('./lib/'))
addpath(genpath('../thirdparty/'))
dataFile = '../data/rooms/office_2p1m_maleSpeech_headMove.mat';
atfFile = '../data/anechoic/glassesOnKemarAnechoicIRs.mat';

%% dereverberation setup
dataStruct = load(dataFile);
arraySigs = dataStruct.arraySigs;
numMics = size(arraySigs,2);

% STFT parameters
blockLenSmp = 2048;
hopSizeSmp = 256;
fftLen = blockLenSmp;
fs = 48000;

% dereverberation parameters
wpeParams.blockLenCovMs = 1500;
wpeParams.adaptiveCovInit = false;
wpeParams.predictionDelayMs = 20; %20; 
wpeParams.regulWeight = 1e-4;

fAx = linspace(0, fs/2, fftLen/2+1)';
wpeParams.predOrders = zeros(fftLen/2+1,1);
wpeParams.predOrders(fAx < 500) = 32;
wpeParams.predOrders(fAx < 1000 & fAx >= 500) = 28;
wpeParams.predOrders(fAx < 8000 & fAx >= 1000) = 24;
wpeParams.predOrders(fAx >= 8000) = 12;

wpeParams.numIterations = 2;
wpeParams.covSmoothingLen = 1;
wpeParams.blockLenStft = blockLenSmp; 
wpeParams.win = sqrt(hann(wpeParams.blockLenStft, 'periodic'));
wpeParams.hopSizeStft = hopSizeSmp;
wpeParams.fftLen = fftLen; 
wpeParams.fs = fs;
wpeParams.numBlocksCovEst = round((wpeParams.blockLenCovMs/1000*fs - (wpeParams.blockLenStft-wpeParams.hopSizeStft)) / wpeParams.hopSizeStft); %292;%floor(3*731/4); % % 344 532 731
wpeParams.overlapFactCov = 0.5;
wpeParams.predDelay = round(wpeParams.predictionDelayMs/1000 / (wpeParams.hopSizeStft / wpeParams.fs)); 


%% dereverberate array signals using WPE
% dereverberate array signals
fprintf('starting multichannel dereverberation...')
arraySigDerev = adaptiveGwpeDerev(arraySigs, wpeParams);
fprintf(' done\n')

noseMicIdx = 4;
fprintf('starting single-channel dereverberation for baseline...')
arraySigDerevBaseline = adaptiveGwpeDerev(arraySigs(:,noseMicIdx), wpeParams);
fprintf(' done\n')

% soundsc(arraySigs(1:5*fs,noseMicIdx),fs)
% soundsc(arraySigDerev(1:5*fs,noseMicIdx),fs)
% soundsc(arraySigDerevBaseline(1:5*fs),fs)

%% apply beamformer
atfDataStruct = load(atfFile);
atfIrsHor = atfDataStruct.irs; % numSamples x numMics x numAtfDirs

% in the struct, we saved rotation angles, if we want to convert that to a
% "measurement grid" we need to flip the sign
atfGridHorAziRad = -atfDataStruct.rotationAnglesDeg*pi/180;

numAtfDirs = size(atfIrsHor,3);
numMics = size(atfIrsHor,2);
numSamplesAtf = size(atfIrsHor, 1);

bfParams.fs = fs;
bfParams.musicFreqRange = [100,1000];
bfParams.bfType = 'MPDR';
bfParams.doaObservationLenMs = 400;
bfParams.regulConst = 1e-2;
bfParams.bfBlockLen = 2048;
bfParams.bfHopSize = bfParams.bfBlockLen/2; 
atfsForBeamformer = fft(atfIrsHor,bfParams.bfBlockLen);
bfParams.useDoaGroundTruth = false;

fprintf('beamforming...')
[arraySigDerevBfd, medianDoaAziEleRad] = musicAndBeamformTimeDomainAdaptive(arraySigDerev, bfParams.bfBlockLen, bfParams.bfHopSize, bfParams.doaObservationLenMs,...
                                            atfsForBeamformer, [atfGridHorAziRad, pi/2*ones(numAtfDirs,1)], bfParams, [], false);
fprintf(' done\n')

% soundsc(arraySigDerevBfd(1:5*fs),fs)

%% setup RIR estimation
chOrder = 2;
numChs = 2*chOrder+1;
chDef = 'complex';

if contains(dataFile,'office')
    blockLenMs = 500; 
    referenceFile = '../data/rooms/office_2p1m_EMA_reference.mat';
elseif contains(dataFile,'storage')
    blockLenMs = 1500;
    referenceFile = '../data/rooms/storage_2p4m_EMA_reference.mat';
elseif contains(dataFile,'kitchen')
    blockLenMs = 1000;
    referenceFile = '../data/rooms/kitchen_2p8m_EMA_reference.mat';
elseif contains(dataFile,'hallway')
    blockLenMs = 2000;
    referenceFile = '../data/rooms/hallway_2p5m_EMA_reference.mat';
end

blockLenSmp = blockLenMs/1000*fs;
fftLen = blockLenSmp;
fAx = linspace(0,fs/2,fftLen/2+1)';

atfs = fft(atfIrsHor,fftLen);

% calculate array matrix from ATFs
encIrLen = 2048;
regulConstRefEncoder = 1e-2;
[D,E,EFilt] = getChdArrayEncDecAtf(atfs, atfGridHorAziRad, fs, chOrder, chDef, fftLen, encIrLen, regulConstRefEncoder);
DOmni = D(:,:,1);

% estimation params
hopSizeSmp = round(blockLenSmp/4);
filterLenSmp = blockLenSmp/2;
win = sqrt(hann(blockLenSmp, 'periodic'));
forgettingFact = 1;
stopEarly = true;

%% convert tracking data
% remove reference pos from position and rotation
translationXYZCart = dataStruct.translationXYZCart - dataStruct.referencePos;
rotYawRad = dataStruct.rotationYawPitchRollRad(:,1) - dataStruct.referenceRotationAziEleRad(1); 
rotEleRad = dataStruct.rotationYawPitchRollRad(:,2) - dataStruct.referenceRotationAziEleRad(2); 

% convert cartesian to radius and angle
[transAziRad,~,transRadius] = cart2sph(translationXYZCart(:,1), translationXYZCart(:,2), zeros(length(translationXYZCart),1));

%% calculate reference SRIR
refStruct = load(referenceFile);
micTypeEma = 'rigid';
c = 343;
chOrderRef = 10;
regulConstRef = 1e-2;

yChRef = getCH(chOrderRef, refStruct.micsAziDeg * pi/180, chDef);
fftLenRefEma = 2^nextpow2(length(refStruct.IRs));
fAxRefEma = linspace(0,fs/2,fftLenRefEma/2+1)';
kRefEma = 2*pi*fAxRefEma/c;

rtfRefEma = fft(refStruct.IRs,fftLenRefEma);
bnRef = circModalCoeffs(chOrderRef, kRefEma*refStruct.arrayRadius, micTypeEma, chDef);

currRadFiltSingleSided = (conj(bnRef) ./ (conj(bnRef) .* bnRef + regulConstRef));
currRadFiltDoubleSided = [currRadFiltSingleSided; conj(flipud(currRadFiltSingleSided(2:end-1,:)))];

chRtfRefEma = currRadFiltDoubleSided .* (pinv(yChRef) * rtfRefEma.').';
rirRefEma = ifft(chRtfRefEma);

%% calculate regularization based on microphone self noise power
noisePowerInterp = interp1(dataStruct.fAxNoise, dataStruct.noisePower, fAx);
regulRls = 1e-2 .* noisePowerInterp;

%% run estimation with reference signal (informed)
fprintf('estimating RIR (informed method)...')
rirEstInformed = fdChRlsRotTransAtf(dataStruct.arraySigs, dataStruct.sourceSig, fs, chOrder, chDef, D, rotYawRad, ...
                                       transAziRad, transRadius, regulRls, forgettingFact, ...
                                       blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly);
fprintf(' done\n')


%% pre-pad the array sigs with zeros to get a causal RIR estimate
prePadLenMs = 5;
prePadLenSmp = prePadLenMs/1000 * fs;

arraySigsPrePad = [zeros(prePadLenSmp, numMics); dataStruct.arraySigs];

minSigLen = min(length(arraySigsPrePad), length(arraySigDerevBfd));
arraySigsPrePad = arraySigsPrePad(1:minSigLen,:);
pseudoRefSigAligned = arraySigDerevBfd(1:minSigLen);
pseudoRefSigAligned(isnan(pseudoRefSigAligned)) = 0; % remove nans (may happen if beamformer has silent input)

rotYawRadPrePad = [zeros(prePadLenSmp, 1); rotYawRad]; 
rotEleRadPrePad = [zeros(prePadLenSmp, 1); rotEleRad]; 
transAziRadPrePad = [zeros(prePadLenSmp, 1); transAziRad];  
transRadiusPrePad = [zeros(prePadLenSmp, 1); transRadius];  
translationXYZCartPrePad = [zeros(prePadLenSmp, 3); translationXYZCart];  

%% run blind estimation (with pseudo reference)
fprintf('estimating RIR (blind)...')
rirEstComp = fdChRlsRotTransAtf(arraySigsPrePad, pseudoRefSigAligned, fs, chOrder, chDef, D, rotYawRadPrePad, ...
                    transAziRadPrePad, transRadiusPrePad, regulRls, forgettingFact, ...
                    blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly);
fprintf(' done\n')

%% blind estimation with omni signal model
fprintf('estimating RIR (blind omni)...')
rirEstOmni = fdChRlsRotTransAtf(arraySigsPrePad, pseudoRefSigAligned, fs, 0, chDef, ...
                DOmni, zeros(size(pseudoRefSigAligned)), zeros(size(pseudoRefSigAligned)), zeros(size(pseudoRefSigAligned)), regulRls, ...
                forgettingFact, blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly);
fprintf(' done\n')

%% baseline
fprintf('estimating RIR (baseline)...')
rirEstBaseline = fdRls(arraySigsPrePad(:,noseMicIdx), arraySigDerevBaseline(1:minSigLen), ...
                fs, forgettingFact, blockLenSmp, hopSizeSmp, filterLenSmp, win, stopEarly);
fprintf(' done\n')

%% plot
ylims = [-40,0];
xlims = [0,0.2*fs];

figure
subplot(511)
plot(db(abs(rirRefEma(:,1))./max(abs(rirRefEma(:,1)))))
grid on
title('reference (from EMA)')
xlim(xlims)
ylim(ylims)

subplot(512)
plot(db(abs(rirEstInformed(:,1))./max(abs(rirEstInformed(:,1)))))
grid on
title('estimate (informed)')
xlim(xlims)
ylim(ylims)

subplot(513)
plot(db(abs(rirEstComp(:,1))./max(abs(rirEstComp(:,1)))))
grid on
title('estimate (blind)')
xlim(xlims)
ylim(ylims)

subplot(514)
plot(db(abs(rirEstOmni)./max(abs(rirEstOmni))))
grid on
title('estimate (blind, omni)')
xlim(xlims)
ylim(ylims)

subplot(515)
plot(db(abs(rirEstBaseline)./max(abs(rirEstBaseline))))
grid on
title('baseline (blind)')
xlim(xlims)
ylim(ylims)

%% resynthesize RIRs 
chOrderResynth = chOrderRef; % use same order as EMA ref
dirSoundResynthAziRad = 0; % direction of direct sound
targetRirLenSmp = round(1.5*fs);
freqRangeSynth = [100,16000]; 

rirEstInformedResynth = synthIsoDiffRirChd(real(rirEstInformed(:,1)), dirSoundResynthAziRad, fs, chOrderResynth, targetRirLenSmp, freqRangeSynth);
rirEstCompResynth = synthIsoDiffRirChd(real(rirEstComp(:,1)), dirSoundResynthAziRad, fs, chOrderResynth, targetRirLenSmp, freqRangeSynth);
rirEstOmniResynth = synthIsoDiffRirChd(real(rirEstOmni), dirSoundResynthAziRad, fs, chOrderResynth, targetRirLenSmp, freqRangeSynth);
rirEstBaselineResynth = synthIsoDiffRirChd(rirEstBaseline, dirSoundResynthAziRad, fs, chOrderResynth, targetRirLenSmp, freqRangeSynth);

figure
subplot(511)
plot(db(abs(rirRefEma(:,1))./max(abs(rirRefEma(:,1)))))
grid on
title('reference (from EMA)')
xlim(xlims)
ylim(ylims)

subplot(512)
plot(db(abs(rirEstInformedResynth(:,1))./max(abs(rirEstInformedResynth(:,1)))))
grid on
title('resynth. estimate (informed)')
xlim(xlims)
ylim(ylims)

subplot(513)
plot(db(abs(rirEstCompResynth(:,1))./max(abs(rirEstCompResynth(:,1)))))
grid on
title('resynth. estimate (blind)')
xlim(xlims)
ylim(ylims)

subplot(514)
plot(db(abs(rirEstOmniResynth(:,1))./max(abs(rirEstOmniResynth(:,1)))))
grid on
title('resynth. estimate (blind, omni)')
xlim(xlims)
ylim(ylims)

subplot(515)
plot(db(abs(rirEstBaselineResynth(:,1))./max(abs(rirEstBaselineResynth(:,1)))))
grid on
title('resynth. baseline (blind)')
xlim(xlims)
ylim(ylims)
