% plotSpecNoise Plots the frequency-domain spectrum in the current axes
%
% This is a QC function to evaluate the consistency of noise across the
% spectrum. A small window is slid accross the spectrum, measuring noiserms
% (and using a polynomial to reduce baseline). When this is plotted, it
% should be relatively flat except for the larger peaks.

% Created: 4/6/2010 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotSpecNoise(sp);

% Operate on the real-valued spectrum
spec = real(sp.spec);
window = 64;

numvals = sp.N - window;
for idx = 1:numvals
   
    shift(idx) = idx-1;
    noiserange = (1:32) + (idx-1);
    specnoise = real(spec(noiserange));
    P = polyfit(noiserange, specnoise',1);
    baselineFit = P(1).*noiserange + P(2);
    specnoise = specnoise - baselineFit';

    noisevar = mean(conj(specnoise) .* (specnoise)) * 2; 
    noiserms(idx) = sqrt(noisevar);
    
end

plot(shift, noiserms);
title('Sliding windown noiserms with polynomial fit');
%xlabel('should be flat (not bowed) except for peaks and near ~1')
%set(gca, 'yscale', 'log')
set(gca, 'ylim', [0 4]);

