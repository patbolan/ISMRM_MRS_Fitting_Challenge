% ****************************************************************************** 
%             MBS - Minnesota Breast Spectroscopy analysis package
% findMax - Finds a maximum point in the frequency domain
%
%   findMax(sp, limits, bAbsMode). Finds the max of the spectrum in real
%   mode (or abs mode, default) and returns the peak ht and position in ppm

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [inten, posppm] = findMax(sp, limits, bAbsMode)


if nargin < 3
    bAbsMode = 1;
end
    
% Get the indices to find a max over
if nargin < 2
    jdxleft = 1;
    jdxright = sp.N;
elseif isempty(limits)
    jdxleft = 1;
    jdxright = sp.N;
else
    % Check for limits at endpoints
    limleft = max(limits(1:2));
    limright = min(limits(1:2));
    
    % Find the indices to search over
    jdxleft = interp1(sp.freqPPM,1:sp.N,limleft,'nearest');
    jdxright = interp1(sp.freqPPM,1:sp.N,limright,'nearest');
    
    if isnan(jdxleft)
        jdxleft = 1;
    end
    if isnan(jdxright)
        jdxright = sp.N;
    end
    
end

for idx = 1:sp.M
    if bAbsMode
        maxjdx = find(abs(sp.spec(:,idx))==max(abs(sp.spec(jdxleft:jdxright,idx))));
    else
        maxjdx = find(real(sp.spec(:,idx))==max(real(sp.spec(jdxleft:jdxright,idx))));    
    end
    
    % Could have a problem if there are two identical points
    if max(size(maxjdx))>1 
        error('findMax failed - found two points with identical values');    
    end
    
    posppm(idx) = sp.freqPPM(maxjdx);
    
    if bAbsMode
        inten(idx) = abs(sp.spec(maxjdx,idx));
    else
        inten(idx) = real(sp.spec(maxjdx,idx));
    end
    
    
end

inten = double(inten);










