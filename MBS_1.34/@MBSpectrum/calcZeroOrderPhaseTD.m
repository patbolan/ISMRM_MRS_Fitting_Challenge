% calcZeroOrderPhaseTD Estimates the zero-order phase of a fid using several
% different algorithms. 
%
%   phs = calcZeroOrderPhaseTD(sp, 1) - uses the most conventional
%   approach; corrects the phase of the first point of the FID (radians)
%
%   phs = calcZeroOrderPhaseTD(sp, N) - uses the first N points of the
%   FID to estimate the first point. Defaults to N=4.
%
% Standard phasing uses only the first point of the FID. This algorithm is
%   an improvement on that. It takes the first N points, does a linear fit,
%   and then uses the predicted value of the phase of the first FID point.

% The results of this phasing depend on the number of points used. Setting
% N=1 gives the conventional algorithm. The default, 4 points, seems to
%   perform well in many cases.

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function phs = calcZeroOrderPhaseTD(sp, N)

if nargin < 2
    N=4;
end

for idx = 1:sp.M
    if N==1
        phs(idx) = angle(sp.fid(1:N,idx));
    else
        phs2fit = unwrap(angle(sp.fid(1:N,idx)));
        p = polyfit(1:N, phs2fit',1);
        
        % Now y=mx+b, x=1. This projects to first FID point.
        phs(idx) = p(1) + p(2);
        %phs(idx) = p(2);        
        
        % TEMP: plot the phase
        %model = p(2) + p(1)*(1:N);
        %plot(1:N, phs2fit, '-k.', 1:N, model, '-r.');
        
    end
        
end























