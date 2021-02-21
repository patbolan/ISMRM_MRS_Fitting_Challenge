% downsample - downsample the FID by integer factors
%
%   sp = downsample(sp, factor) 
% Used for data (eg Siemens) that were oversampled. 
%   This works by a hard truncation in the spectral domain, which leads to an
%   artifact at the end of the time-domain fid. 
% Alternatively, try truncateFreqDomain(). This version does cause a slight
% frequency shift that is not properly corrected, whereas
% trucateFreqDomain() does not.

% Created: 1/2/2009 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function sp = downsample(sp, factor)

% default to 2x
if nargin < 2
    factor = 2;
end

% Copy pts from the old spec to the new one
newpts = sp.N / factor;
newpts = floor(newpts/2) * 2; % Better be divisible by 2
edge = (sp.N - newpts)/2; % # pts to remove on each edge

% Just keep the middle freq domain points
sp.spec = sp.spec(edge+1:sp.N-edge, :);

% Adjust the spectral width
sp.sw = sp.sw/ factor;












