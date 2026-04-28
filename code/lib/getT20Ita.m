function [t20, freqAxis, rirEndSmp] = getT20Ita(rir, fs, bandsPerOctave, freqRange, fillMissingVals)
% calculates reverberation time as average of T15, T20, T30
arguments
    rir
    fs
    bandsPerOctave
    freqRange
    fillMissingVals = true
end

itaIr = itaAudio;
itaIr.samplingRate = fs;
itaIr.time = rir;
itaTn = ita_roomacoustics(itaIr, 'freqRange', freqRange, 'bandsPerOctave', bandsPerOctave, 'T20', 'T15', 'T10', 'Intersection_Time_Lundeby'); % 'plotLundebyResults'
t20 = itaTn.T20.freqData;
t15 = itaTn.T15.freqData;
t10 = itaTn.T10.freqData;

% replace missing values of T20 with T15 or T10 results
t20(isnan(t20)) = t15(isnan(t20));
t20(isnan(t20)) = t10(isnan(t20));

if fillMissingVals
    while any(isnan(t20), 'all')
        warning('getT20Ita: found nans in RT estimation, interpolating')
        % rt60 = fillmissing(rt60,'linear'); % linear interpolation may result in negative RT60s
        t20 = fillmissing2(t20,"movmedian",3); % replace missing values by the 3x3 median, repeat process if necessary
    end
end

freqAxis = itaTn.T20.freqVector;
intersectionTimeLundeby = itaTn.Intersection_Time_Lundeby.freqData;

rirEndSmp = round(mean(intersectionTimeLundeby * fs, "all"));
rirEndSmp = min(rirEndSmp, size(rir,1));


