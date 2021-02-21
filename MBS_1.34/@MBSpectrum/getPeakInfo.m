% getPeakInfo - calculates information on a peak withing region
%
%  [inten, posppm, corrected_inten, ...
%   fwhm_ppm, posppm_left_hm, posppm_right_hm, area, corrected_area ] = ...
%    getPeakInfo(sp, limits, bAbsMode)
% Calculates several parameters for a the largest "peak" located with the
% region. 
% inten = intensity (height)
% posppom = position of max 
% corrected_inten = intensity, but corrected using a linear baseline
%   between the two limits
% fwhm_ppm = full width, half max, in ppm
% area = area of entire region between limits
% corrected_area = area, but corrected using a linear baseline
%   bAbsMode - 1 for abs mode, 0 for real mode

% Created: 2/6/2006 Patrick Bolan,
% MBS - Minnesota Breast Spectroscopy package
function [inten, posppm, corrected_inten, fwhm_ppm, ...
    area, corrected_area ] = ...
    getPeakInfo(sp, limits, bAbsMode)

if nargin < 3
    bAbsMode = 1;
end

% Get the indices to find a max over
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
    
    % Check th ends
    if (jdxleft < 1) || ~isfinite(jdxleft)
        jdxleft = 1;
    end
    if (jdxright>sp.N) || ~isfinite(jdxright)
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

    % Now calculate FWHM. I originally used the interp1 function, which is
    % fast, but the reverse of a spectrum is not a function (ie, a given
    % amplitude does not uniquely define a frequency), so it is
    % error-prone. Instead I'm performing a tedious walk-down from the max
    % point.
    if bAbsMode
        inten_left(idx) = abs(sp.spec(jdxleft,idx));
        inten_right(idx) = abs(sp.spec(jdxright,idx));
        spec_frag = abs(sp.spec(:,idx));
    else
        inten_left(idx) = real(sp.spec(jdxleft,idx));
        inten_right(idx) = real(sp.spec(jdxright,idx));
        spec_frag = real(sp.spec(:,idx));
    end
    corrected_inten(idx) = inten(idx) - (inten_left + inten_right)/2;

    % This is the target half-height in each direction
    halfht_left = inten(idx) - (inten(idx) - inten_left(idx))/2;
    halfht_right = inten(idx) - (inten(idx) - inten_right(idx))/2;

    % Walk left and linearly interpolate to find the half-ht
    curjdx = maxjdx - 1;
    posppm_left_hm = sp.freqPPM(curjdx);
    
    while (curjdx >= jdxleft)
        if(spec_frag(curjdx) <= halfht_left)
            % We've gone below the half-ht. Interp to find position
            posppm_left_hm(idx) = interp1(spec_frag(curjdx:curjdx+1), ...
                sp.freqPPM(curjdx:curjdx+1),halfht_left, 'linear');
            break;
        end
        curjdx = curjdx - 1;
    end

    % Start at the max point
    curjdx = maxjdx + 1;
    posppm_right_hm(idx) = sp.freqPPM(curjdx);
    while (curjdx <= jdxright)
        if(spec_frag(curjdx) <= halfht_right)
            % We've gone below the half-ht. Interp to find position
            posppm_right_hm(idx) = interp1(spec_frag(curjdx:-1:curjdx-1), ...
                sp.freqPPM(curjdx:-1:curjdx-1),halfht_right, 'linear');
            break;
        end
        curjdx = curjdx + 1;
    end
    %posppm_left_hm = interp1( spec_fragment, freq_fragment, halfht_left, 'linear');
    %posppm_right_hm = interp1( spec_fragment, freq_fragment, halfht_right, 'linear');

    fwhm_ppm(idx) = posppm_left_hm(idx) - posppm_right_hm(idx);
    
    % Integrals
    % Units are per-point, not per-ppm
    area(idx) = trapz(spec_frag(jdxleft:jdxright));
    corrected_area(idx) = area(idx) - ...
        (inten_left(idx) + inten_right(idx)) * (jdxright - jdxleft) / 2;
    
    % Correct areas so that they are in true units, ppm*ht
    ppm_per_point = sp.freqPPM(1) - sp.freqPPM(2);
    area(idx) = area(idx) * ppm_per_point;
    corrected_area(idx) = corrected_area(idx) * ppm_per_point;  
end












