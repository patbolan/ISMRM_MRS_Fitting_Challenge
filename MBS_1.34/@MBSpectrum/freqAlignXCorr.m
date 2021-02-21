% freqAlignXCorr - Align spectra array using cross-correlation
%
%  [sp, shiftppm] = freqAlignXCorr(sp, range, resolutionHz)
%   Shifts the spectra. Finds the max point between limits, and sets that peak 
%   to the reference ppm specified.
% To calculate the x-correlation, it will compare each spectrum to the mean
% spectrum, and determine the max of the xcorrelation. 

% Created: 7/15/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [sp, shiftppm] = freqAlignXCorr(sp, range, resolutionHz)

numspec = sp.M;
sfrq = sp.frq;
swppm = sp.swppm;
pts = sp.N;
deltaFppm = sp.swppm / (sp.N-1);
deltaFHz = deltaFppm * sp.frq;



% Pre-calculate how much ZF is necessary
zf = ceil(deltaFHz/resolutionHz);
zf = max(zf,1);
newDeltaFHz = deltaFHz/zf;
disp(sprintf('Native FD res %.2f Hz, requested %.2f Hz, need zf %.0f gives %.2f Hz',...
    deltaFHz, resolutionHz, zf, newDeltaFHz));

% Pre-calculate Range
maxlags = ceil(range/newDeltaFHz);
disp(sprintf('To cover %.2f Hz, need %.0f points', range, maxlags));

% Zero-fill first
shiftppm = zeros(numspec,1);
sp_sm = zeroFill(sp, pts*zf);
spec = abs(sp_sm.spec);
sp_mean = mean(sp_sm);
specref = abs(sp_mean.spec);


% Loop over each.
for idx = 1:sp.M
    %s1 = spec(:,1);
    s1 = specref;
    s2 = spec(:,idx);
    c = xcorr(s1, s2, maxlags);
    
    % Find the max point
    peakpt = find(c==max(c));
    offsetIdx(idx) = peakpt - (maxlags+1);
    shiftHz(idx) = offsetIdx(idx) * newDeltaFHz;
   
    %fprintf('%d: %.3f Hz\n', idx, shiftHz(idx));
    
%     % DEBUG
%     figure(202)
%     xcrng = -maxlags:1:maxlags;
%     plot(xcrng,c);
%     
%     freq = sp_sm.freq;
%     sfrq = sp_sm.frq;
%     freq = freq .* sfrq;
%     figure(201)
%     plot(freq, s1, '-k',...
%         freq, s2, '-b');
%     set(gca,'xlim',[0 1000]);
%     legend(sprintf('%.0f', idx-1), sprintf('%.0f', idx));
%     %dummy = input(sprintf('Shift %.0f pts, %.2f Hz', offsetIdx(idx), shiftHz(idx)));
    
end


% Check that the range was sufficient.
if (max(abs(offsetIdx(:))) == maxlags)
    % The max was at the end. Probbably not OK.
    % be more robust - allow 1 or two to be this large.
    disp(sprintf('freqAlignXCorr: Likely insufficient Freq range'))
end

% The shifts shiftHz are differential. Convert them to absolute
% shiftHzAbs(1) = shiftHz(1);
% for idx = 2:numspec
%     shiftHzAbs(idx) = shiftHz(idx) + shiftHzAbs(idx-1);    
% end
shiftHzAbs = -shiftHz;
%shiftHzAbs = fliplr(shiftHz);

% Now apply the shifts, in ppm
shiftppm = shiftHzAbs ./ sfrq;

sp = shiftFreq(sp, shiftppm');











