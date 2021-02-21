% hsvdRemoveBroadPeaks - Removes any peaks with linewidths greater than
% specified
%
%   [residual, peaks, spfit] = hsvdRemoveOutside(sp, fwhm, maxpts)
%       sp - input MBSpectrum
%       fwhm - the cuttoff full-width-half-max, in Hz. Anthing larger is
%       removed
%       maxpts - max number of points to use when calculating SVD. Defaults
%       to using all, but the calculation scales with N^2
%       residual - an MBSpectrum of the residual
%       peaks - an array of lorentzian peak values. See simulateLorentzian
%           for structure
%       spfit - an NxM MBSpectrum of all M singular values
%

% Created: 8/11/2010 Patrick Bolan and Timo Liimatainen
% MBS - Minnesota Breast Spectroscopy package
function [residual, peaks, spfit] = hsvdRemoveBroadPeaks(sp, fwhm, maxpts)

% Keep the original sp for final calculation, but do SVD on sp, which is
% truncated in the time domain. This works pretty well since the signal
% information is weighted to the beginning of the FID (especially broad
% peaks)
sporig = sp;
if nargin>2
    sp = sporig.resizeTimeDomain(maxpts);
end

% First, normal HSVD
[~, peaks, ~] = hsvd(sp);

badIdx = [];
for idx=1:size(peaks,2);
    
    if (peaks(idx).lambda > fwhm*pi)
        badIdx = [badIdx idx]; 
    end
end

peaksToRemove = peaks(badIdx);

spfit = simulateLorentzian(peaksToRemove, sporig.N, sporig.sw, sporig.frq);
residual = sporig - spfit.sum;






