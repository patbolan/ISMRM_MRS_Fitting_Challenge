% gaussianFilter - Gaussian line-broadening
%
% sp = gaussianFilter(sp, gaussLW) - Multiplies the FID by an squared 
% exponential function. 
%   gaussLW is the linebroadening characteristic in Hz. 
% gaussLW is equivalent to the full-width-half-max (FWHM). For a pure 
% gaussian line, adding this filter will increase the FWHM by gaussLW.
% The gauss linewidth is not directly related to R2star
%
%   exp(-(gamma*t)^2), where gamma = gaussLW * pi /(2 sqrt(ln(2)))

% Note also the relationships:
%   lorentzLW = FHWM
%   lorentzLW = 1/(pi*t2)
%   lambda = 1/t2 = R2star


% Created: 2060701 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function sp = gaussianFilter(sp, gaussLW)

gamma = gaussLW .* (pi/(2*sqrt(log(2))));
sp.fid = sp.fid .* (exp(-(sp.time * gamma).^2) * ones(1,sp.M));











