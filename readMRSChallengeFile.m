function [sp, spW] = readMRSChallengeFile(fname)


% Cols are real img real imag
fp = fopen(fname, 'r');
C = textscan(fp, '%f%f%f%f');
fclose(fp);


fidM = C{1} + 1j*C{2};
fidW = C{3} + 1j*C{4};

% Note - the FIDs from gosia have the first point multiplied by 0.5
% already! Fix this



sp = MBSpectrum;
sp.frq = 123.2;
sp.sw = 4000;
sp.ppmref = 0;
sp.params.EffectiveEchoTime.value{1} = num2str(0.030);
sp.params.EffectiveEchoTime.unit = 'ms';


sp.fid = fidM;

spW = sp;
spW.fid = fidW;
fprintf('read %d points.\n', sp.N);