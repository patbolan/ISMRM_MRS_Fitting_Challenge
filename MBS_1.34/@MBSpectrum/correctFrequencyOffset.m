% Analogous to correctDC, this looks for an offset in the frequency domain
% and subtacts it off.
% This is similar to scaling the the first point of the FID, and is
% described in Zhu et al, JMR A 105, 1993

% Created: 2/28/2011 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [sp_mod, correctionFD] = correctFrequencyOffset(sp)

% Some fixed parameters
prefilterLB = 8;
halfrange = 10;

sp_mod = [];

correctionFD = zeros(1, sp.M);
for idx=1:sp.M
    
    spone = sp.extract(idx);
    
    [tmp, correctionFD] = correctFrequencyOffsetSingle(spone);
    if isempty(sp_mod)
        sp_mod = tmp;
    else
        sp_mod = sp_mod.append(tmp);
    end
    
end



function [spmod, correctionFD] = correctFrequencyOffsetSingle(sp)

spec = sp.spec;

rangemin = 1:20;

% Take the complex average over this range
correctionFD = mean(spec(rangemin));
spec = spec - correctionFD;

% Now assign the output
spmod = sp;
spmod.spec = spec;







