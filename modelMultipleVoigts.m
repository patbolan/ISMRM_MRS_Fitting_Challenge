% This function will model a spectrum with any number of Voigt resonances,
% and return the modeled spectrum, as well as the real-valued trace
% interpolated onto the frequency values provided.
%
% beta - the Voigt parameters. Each peak has 5, and they should be in order
% (ie if you have 3 peaks there should be 15 params). The parameters will
% be interpreted as:
%   beta(1) - amplitude, arbitrary units
%   beta(2) - center frequency, ppm
%   beta(3) - phase, degrees
%   beta(4) - lorentzian linewidth (FWHM), Hz
%   beta(5) - guassian FWHM, Hz
%
% sp_ref is a reference spectum, defining time and frequency axes
% freqvals are the values for sampling the resultant realvalued spectra
%   onto. This is used for fitting by lsqcurvefit - you can just select those
%   frequency values you want to fit over. Don't need to be continuous.
% sp_model is the resultant fit spectrum.
%
% Here's how to test with a simple water & lipid
% F = modelMultipleVoigts([1, 4.7, 0, 5, 5, 2, 1.3, 0, 20, 0], 0:.01:6, sp_ref);
% 20110301 PJB
function [simReSpec, sp_model, sp_each] = modelMultipleVoigts(beta, freqvals, sp_ref)

% The procedure is to calculate the FID analytically, FFT to frequency
% domain, and then numerically interpolate onto the specified frequencies.

% beta has 5 values for each peak 
numpeaks = floor(size(beta,2) / 5);
npeach = 5;

% Extract data from the reference spectrum
timeax = sp_ref.time;
freqax = sp_ref.freqPPM;
sfrq = sp_ref.frq;
ppmref = sp_ref.ppmref;

% First, simulate the time domain data
pts = max(size(timeax));
modelFID = zeros(size(timeax,1), 1);
sp_each = [];
for idx = 1:numpeaks
    
    amp = beta(1+(idx-1)*npeach);
    frqppm = beta(2+(idx-1)*npeach);
    phs = beta(3+(idx-1)*npeach);
    llw = beta(4+(idx-1)*npeach);
    glw = beta(5+(idx-1)*npeach);
    
    %   Lorentz:  exp(-lambda*t);     lambda = FWHM * pi
    %   Gauss:    exp(-(gamma*t)^2);  gamma = FWHM * pi / (2 sqrt(ln(2)) )
    lambda = llw .* pi;
    gamma = glw .* (pi/(2*sqrt(log(2))));
    frq = -(frqppm - ppmref) * sfrq;

    modelFID = amp .* exp( ...
        (1i * 2 * pi *  frq) .* timeax + ...
        (1i * phs / 180 * pi) + ...
        (-lambda .* timeax) + ...
        (-(gamma.*timeax).^2)); 
   
    
    % Create the model spectrum
    sp_tmp = sp_ref;
    sp_tmp.fid = modelFID;
    
    % Append to array
    if isempty(sp_each)
        sp_each = sp_tmp;
    else
        sp_each = sp_each.append(sp_tmp);
    end
    
end

% Need one result for the lsqnonlin 
sp_model = sp_each.sum;

% FFT, being smart about the first point and scaling
tmp = sp_model.fid;
tmp(1) = tmp(1) * 0.5;
spec_c = fftshift(fft(tmp), 1);
spec_c = spec_c ./ sqrt(pts);
spec_r = real(spec_c);

% Plus, a vertical offset is the last beta
%spec_r = spec_r + beta(6+(numpeaks-1)*npeach);


% Now interpolate these data onto the frequency values requested.
simReSpec = interp1(freqax, spec_r, freqvals, 'pchip');


% DEBUG - plot
% figure(101)
% subplot(3,1,1)
% plot(timeax, real(modelFID));
% ylabel('SimFid');
% subplot(3,1,2)
% plot(freqax, spec_r); 
% set(gca, 'xdir', 'reverse');
% ylabel('SimSpec');
% subplot(3,1,3)
% plot(freqvals, simSpec);
% ylabel('SimSpec, sampled');
% set(gca, 'xdir', 'reverse');
