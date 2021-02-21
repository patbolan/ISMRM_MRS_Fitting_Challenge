% This function expects a direcotry called "basisset_text" on the current path

function [spBasis, labels] = readMRSChallengeBasisSet()

dirname = 'basisset_text';

% Hardwire stuff, quick and dirty
labels = {'Ace', 'Ala', 'Asc', 'Asp', 'Cr', 'GABA', 'Glc', ...
    'Gln', 'Glu', 'Gly', 'GPC', 'GSH', 'Ins', 'Lac', 'Mac', ...
    'NAAG', 'NAA', 'PCho', 'PCr', 'PE', 'sIns', 'Tau' };
labels = {'Mac', 'Ala', 'Asc', 'Asp', 'Cr', 'GABA', 'Glc', ...
    'Gln', 'Glu', 'Gly', 'GPC', 'GSH', 'Ins', 'Lac', ...
    'NAAG', 'NAA', 'PCho', 'PCr', 'PE', 'sIns', 'Tau' };
spBasis = MBSpectrum;
spBasis.frq = 123.2;
spBasis.sw = 4000;
spBasis.ppmref = 0;
spBasis.params.EffectiveEchoTime.value{1} = num2str(0.030);
spBasis.params.EffectiveEchoTime.unit = 'ms';

Np = 2048;
fid = zeros(Np, size(labels,2));
for idx=1:size(labels,2)
    fname = fullfile(dirname, sprintf('%s.txt', labels{idx}));
    fprintf('reading %s\n', fname);
    fp = fopen(fname, 'r');
    C = textscan(fp, '%f%f');
    fclose(fp);
    
    fid(:, idx) = C{1} + 1j*C{2};    
end
spBasis.fid = fid;

fprintf('complete Basis set.\n');






