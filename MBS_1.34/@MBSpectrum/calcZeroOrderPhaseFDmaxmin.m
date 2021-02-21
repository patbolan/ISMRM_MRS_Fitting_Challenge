% calcZeroOrderPhaseFDmaxmin Estimates the zero-order phase of a fid by
% maximizing the minimum values in the frequency domain.
%
% phs = calcZeroOrderPhaseFDmaxmin(sp, limits) - Estimates
%   the zero order phase. The limits (in ppm, optional) can be used to define a 
%   section of the frequency domain
%
% Finds the zeroth order phase of the spectrum, and returns the angle in 
%   radians.
% A Frequency Domain algorithm - finds the zero-order phase
%   which maximizes the minimum of the real spectrum. First, it does some
%   linebroadening on the peak because it is susceptable to noise.
% This one is particularly helpful when working with spectra with lipid
% resonances in it. Other methods get confused by the asymmetric
% resonances, and can find poor optima.
% Works only on a single fid at a time. Works in three passes. 

% Created: 5/23/2003 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function phs = calcZeroOrderPhaseFDmaxmin(sp)

% This is susecptabile to noise; smooth a little first
sp = lorentzianFilter(sp, 2);

% Lookup frequency-domain limits
if nargin < 2
   limits = [max(sp.freqPPM(:)), min(sp.freqPPM(:))];
end

phs = zeros(1, sp.M);

for jdx = 1:sp.M
    spSingle = extract(sp,jdx);

    % Pass #1: find the global max. 
    %rp = 0;
    %incr  = 2*pi/20;
    %rp = findOptimalPhase(sp, limits, rp, incr);

    % PJB 20100903 Modification: for first pass, use simpler (and faster)
    %   time-domain estimation
    rp = -1 * calcZeroOrderPhaseTD(spSingle, 4);
    
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
    minspec(idx) = min(spec(:));
end
rp = rps(find(minspec==max(minspec)));
rp = rp(1);
%plot(rps, minspec);











