% ****************************************************************************** 
%             MBS - Minnesota Breast Spectroscopy analysis package
% findMin - Finds a maximum point in the frequency domain
%
%   findMin(sp, limits, bAbsMode). Finds the max of the spectrum in real
%   mode (or abs mode, default) and returns the peak ht and position in ppm

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [inten, posppm] = findMin(sp, limits, bAbsMode)


if nargin < 3
    bAbsMode = 1;
end
    
% Get the indices to find a min over
if nargin < 2
    jdxleft = 1;
    jdxright = sp.N;
else
    % Check for limits at endpoints
    limleft = max(limits(1:2));
    limright = min(limits(1:2));
    
    % Find the indices to search over
    jdxleft = interp1(sp.freqPPM,1:sp.N,limleft,'nearest');
    jdxright = interp1(sp.freqPPM,1:sp.N,limright,'nearest');
end

for idx = 1:sp.M
    if bAbsMode
        minjdx = find(abs(sp.spec(:,idx))==min(abs(sp.spec(jdxleft:jdxright,idx))));
    else
        minjdx = find(real(sp.spec(:,idx))==min(real(sp.spec(jdxleft:jdxright,idx))));    
    end
    
    % Could have a problem if there are two identical points
    if min(size(minjdx))>1 
        error('findMin failed - found two points with identical values');    
    end
    
    posppm(idx) = sp.freqPPM(minjdx);
    
    if bAbsMode
        inten(idx) = abs(sp.spec(minjdx,idx));
    else
        inten(idx) = real(sp.spec(minjdx,idx));
    end
    
    
end

inten = double(inten);










