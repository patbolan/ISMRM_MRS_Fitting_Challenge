% getSnrFd - Returns the snr, noiserms, and maxsig in the Frequency domain
% 	for the summed Spectra.
%
%   [snr_amp, snr_pwr, noiserms, maxsig] = getSnrFd(sp)
%   Noise is always for complex data = 2x real data
%   This function returns a simple ratio - you can convert it to dB as follows:
%       snr dB  = 10*log10(snr_amp) 
%               = 20*log10(snr_pwr)
%
%  NOTE: This assumes that there are no signals on the left part of the
%  frequency spectrum. This is not valid if there are artifacts or aliased
%  peaks there (ie, EPSI). Also assumes white (frequency-independent)
%  noise, which is not the case with Siemens data unless the oversampling
%  is removed. 

% Created: 1/18/2003 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [snr_amp, snr_pwr, noiserms, maxsig] = getSnrFd(sp)

spave = mean(sp);
spec = spave.spec;

noiserange = 1:1/16 * sp.N;

% I calculate the noise with real data. I could 
%   do real and imaginary separately and add them 
%   together, but they should be equal, so I'll 
%   do the real only and double it.
specnoise = real(spec(noiserange));
P = polyfit(noiserange, specnoise',1);
baselineFit = P(1).*noiserange + P(2);
specnoise = specnoise - baselineFit';
    
noisevar = mean(conj(specnoise) .* specnoise) * 2; 
noiserms = sqrt(noisevar);
maxsig = abs(max(spec));
snr_amp = maxsig / (noiserms); 

% Also snr_pwr
snr_pwr = (1/sp.N)*sum(spec.*conj(spec))/(noisevar);
