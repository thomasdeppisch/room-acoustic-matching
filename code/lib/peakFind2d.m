function peaks_idx = peakFind2d(pmap, grid_xyz, nPeaks)
% function [peak_idx, peak_dirs] = peakFind2d(pmap, grid_dirs, grid_xyz, nPeaks)

kappa  = 50;
P_minus_peak = pmap;
%peaks_dirs = zeros(nPeaks, 2);
peaks_idx = zeros(nPeaks, 1);
for k = 1:nPeaks
    [~, peak_idx] = max(P_minus_peak);
%    peaks_dirs(k,:) = grid_dirs(peak_idx,:);
    peaks_idx(k) = peak_idx;
    if (k~=nPeaks)
        VM_mean = grid_xyz(peak_idx,:); % orientation of VM distribution
        VM_mask = kappa/(2*pi*exp(kappa)-exp(-kappa)) * exp(kappa*grid_xyz*VM_mean'); % VM distribution
        VM_mask = 1./(0.00001+VM_mask); % inverse VM distribution
        P_minus_peak = P_minus_peak.*VM_mask;
    end
end
end
