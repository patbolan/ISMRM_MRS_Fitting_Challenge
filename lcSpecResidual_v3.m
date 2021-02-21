% v3 is v2, but switching to Voigt filtering
% v2 is tailored just for the iftting contest. 
% Includes amp, phase, r2s, delf for each 

function [residual, spSyn] = lcSpecResidual_v3(theta, spdata, spbasis, rangePPM)

if ~exist('rangePPM', 'var')
    freqAx = spdata.freqPPM;
    rangePPM = [min(freqAx(:)) max(freqAx(:))];
end
    
J = spbasis.M;
nGlobals = 4;
nParamsPerElement = 5;

% Synthesize
for jdx=1:J
    joff = nGlobals + nParamsPerElement*(jdx-1);
    
    % get per-basis parameters
    amp = theta(joff+1);
    phs = theta(joff+2);
    llw = theta(joff+3);
    glw = theta(joff+4);
    delf = theta(joff+5);
       
    % Units
    freqppm = delf/spdata.frq;
    phsrad = phs .* pi/180;
    
    spone = spbasis.extract(jdx);
    spone = spone.times(amp);
    spone = spone.phaseShift(phsrad);
    spone = spone.lorentzianFilter(llw);
    spone = spone.gaussianFilter(glw);
    spone = spone.shiftFreq(freqppm);
        
    if jdx==1
        spSyn = spone;
    else
        spSyn = spSyn + spone;
    end
end

% Global corrections
phsGlb = theta(1);
llwGlb = theta(2);
glwGlb = theta(3);
delfGlbHz  = theta(4);
delfGlbPPM  = delfGlbHz / spdata.frq;

spSyn = spSyn.shiftFreq(delfGlbPPM);
spSyn = spSyn.phaseShift(-phsGlb);
spSyn = spSyn.lorentzianFilter(llwGlb);
spSyn = spSyn.gaussianFilter(glwGlb);


% calculate residual 
spdiff = spdata - spSyn;
%residual = [real(spdiff.spec) imag(spdiff.spec)];

% Frequency-specific
subspec = spdiff.getSubSpec(rangePPM);
residual = [real(subspec) imag(subspec)];


% Very slow debug! 
figure(19)
freqax = spdiff.freqPPM;
plot(freqax, real(spdata.spec), 'k', ...
    freqax, real(spSyn.spec), 'b', ...
    freqax, real(spdiff.spec), 'r');
set(gca, 'xdir', 'reverse')
%legend('data', 'fit', 'residual')
%zoom on
set(gca, 'xlim', [rangePPM(1)-1 rangePPM(2)+1]);
% 
% figure(20)
% timeax = spdiff.time;
% plot(timeax, real(spdata.fid), 'k', ...
%     timeax, real(spSyn.fid), 'b', ...
%     timeax, real(spdiff.fid), 'r');

%legend('data', 'fit', 'residual')
%zoom on
pause(0.01);

return;




