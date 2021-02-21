% getSnrTd - Returns the snr, noiserms, and maxsig in the Time domain
% 	for the summed Spectra.
%
%   [snr_amp, snr_pwr, noiserms, maxsig] = getSnrTd(sp)
%   Noise is always for complex data = 2x real data
%   This function returns a simple ratio - you can convert it to dB as follows:
%       snr dB  = 10*log10(snr_amp) 
%               = 20*log10(snr_pwr)
%
%  NOTE: This assumes that there is no coherent signal left at the end of
%  the spectrum. This is invalid for phantoms, long T2*, short
%  acquisitions, etc

% 20160623 PJB Adding polynomial fitting just as in the getSnrFD. this is
% important for phantoms.

% Created: 1/18/2003 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [snr_amp, snr_pwr, noiserms, maxsig] = getSnrTd(sp)

% Best to do a correctDC to get the noise at the end averaging 0
sp = correctDC(sp);

sp = mean(sp);
noiserange = 7/8 * sp.N: sp.N;

fidnoise = real(sp.fid(noiserange));

% Poly fitting to remove baseline
P = polyfit(noiserange, fidnoise',1);
baselineFit = P(1).*noiserange + P(2); % 1st order polynomial fitting
fidnoise = fidnoise - baselineFit';
noisevar = mean(conj(fidnoise) .* fidnoise) * 2;
%noisevar = mean((conj(fid(noiserange)) .* fid(noiserange))); % Complex noise

noiserms = sqrt(noisevar);
maxsig = abs(max(sp.fid)); % Max complex signal
snr_amp = maxsig / (noiserms); 


% I should be calculating the SNR using a power measure, not just an amplitude 
%   ratio. Assume the noise value is correct. Then use the power for the 
%   signal power. 
snr_pwr = (1/sp.N)*sum(sp.fid.*conj(sp.fid))/(noisevar);
