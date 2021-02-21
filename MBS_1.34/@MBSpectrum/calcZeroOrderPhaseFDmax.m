% calcZeroOrderPhaseFDmax Estimates the zero-order phase of a fid by
% maximizing the values in the frequency domain.
%
% phs = calcZeroOrderPhaseFDmax(sp, limits) - Estimates
%   the zero order phase. The limits (in ppm) can be used to define a 
%   section of the frequency domain
%
% A Frequency Domain algorithm - finds the zero-order phase
%   which maximizes the sum of the real spectrum. 
% This method is pretty much unaffected by smoothing.
% Works only on a single fid at a time. Uses degrees internally, but 
%   inputs/outputs are all radians 
% Algorithm: A crude optimization technique - maximizes the real part of
%   the FD spectrum. Works in three passes of increasing resolution
% Note that this is slow, now that the spec is calculated each time rather
% than stored directly as in previous versions. To speed this up, extract
% the spectral values and do the phasing explicitly in this function.

% Created: 5/23/2003 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function phs = calcZeroOrderPhaseFDmax(sp, limits)

% Lookup limits
if nargin < 2
   limits = [max(sp.freqPPM(:)), min(sp.freqPPM(:))];
end

phs = zeros(1, sp.M);

for jdx = 1:sp.M
    spSingle = extract(sp,jdx);

    % Pass #1: find the global max. 
    rp = 0;
    incr  = 2*pi/20;
    rp = findOptimalPhase(sp, limits, rp, incr);

    % Pass #2: Find the phase to 1-degree accuracy. 
    clear realsumsq;
    incr = 1 * (pi/180);
    rp = findOptimalPhase(sp, limits, rp, incr);
    
    % Pass #3: 0.1 degree accuracy.
    incr = 0.1 * (pi/180);
    rp = findOptimalPhase(sp, limits, rp, incr);
    
    % Pass #4: 0.01 degree accuracy.
    incr = 0.01 * (pi/180);
    rp = findOptimalPhase(sp, limits, rp, incr);

    
    % Return the negative; this is the amount needed to correct the current
    % zero-order phase
    phs(jdx) = -rp(1);   
end


% This function generates a new set of phases around rp and tests,
% returning the optimal value
function rp = findOptimalPhase(sp, limits, rp, incr)
rps = (rp-10*incr:incr:rp+10*incr);
for idx = 1:size(rps,2)
    temp = truncateFreqDomain(phaseShift(sp,rps(idx)), limits);
    spec = real(temp.spec);
    realsumsq(idx) = sum(spec(:));
end
rp = rps(find(realsumsq==max(realsumsq)));
rp = rp(1);
%plot(rps, minspec);






