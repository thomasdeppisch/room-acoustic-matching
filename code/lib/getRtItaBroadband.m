function [t20,rirEndSmp] = getRtItaBroadband(rir, fs, freqRange)
% calculates reverberation time as average of T15, T20, T30

itaIr = itaAudio;
itaIr.samplingRate = fs;
itaIr.time = rir;
itaTn = ita_roomacoustics(itaIr, 'freqRange', freqRange, 'broadbandAnalysis', true, 'T20', 'T15', 'T10', 'Intersection_Time_Lundeby'); % 'plotLundebyResults'
t20 = itaTn.T20.freqData;
t15 = itaTn.T15.freqData;
t10 = itaTn.T10.freqData;

% replace missing values of T20 with T15 or T10 results
t20(isnan(t20)) = t15(isnan(t20));
t20(isnan(t20)) = t10(isnan(t20));

freqAxis = itaTn.T20.freqVector;
intersectionTimeLundeby = itaTn.Intersection_Time_Lundeby.freqData;

rirEndSmp = round(mean(intersectionTimeLundeby * fs, "all"));
rirEndSmp = min(rirEndSmp, size(rir,1));


