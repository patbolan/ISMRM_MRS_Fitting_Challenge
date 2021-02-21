% lorentzianFilter - Lorentzian line-broadening
%
% sp = lorentzianFilter(sp, lorentzLW) - Multiplies the FID by an exponential. 
%   lorentzLW is the linebroadening characteristic in Hz, lb. 
% lorentzLW is equivalent to the full-width-half-max (FWHM). For a pure 
% Lorentzian line, adding this filter will increase the FWHM by lorentzLW
% Note that lorentzLW = FHWM = 1/(pi*t2)

% 20160701 PJB: I noticed an error in this code. It acutally doubled the
% filterwidth! In correct notation, the filter is 
%   exp(-lambda*t), where lambda = lorentzLW * pi 
% Note also the relationships:
%   lorentzLW = FHWM
%   lorentzLW = 1/(pi*t2)
%   lambda = 1/t2 = R2star


% Created: 8/7/2002 Patrick Bolan, with revisions by Curt Corum
% MBS - Minnesota Breast Spectroscopy package
function sp = lorentzianFilter(sp, lorentzLW)

lambda = lorentzLW * pi;
sp.fid = sp.fid .* (exp(-sp.time * lambda) * ones(1,sp.M));











