% go_fit_v4 -
% This is my demo fit code for publishing

% I did these cases one at a time you can loop if you prefer
dsNum = 1;
fprintf('Processing case #%d.\n');

% Options, can vary for each case
bPrePhaseWater = false;
bHSVDWaterRemoval = false;
bIncludeLipid = false;

%%
outdir = fullfile(pwd, 'output');
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
fnameResults = fullfile(outdir, sprintf('results_%02d.txt', dsNum)); 
fnameWaterFit = fullfile(outdir, sprintf('waterfit_%02d.jpg', dsNum)); 
fnameMetaboliteFit = fullfile(outdir, sprintf('metabolitefit_%02d.jpg', dsNum)); 


dirname = fullfile(pwd, 'datasets_text');
fname = fullfile(dirname, sprintf('dataset%d.txt', dsNum));
[sp, spW] = readMRSChallengeFile(fname);

% The metabolite spectra have their first point is already scaled
% down by 0.5. Water does not. MM says its a LCModel thing
sp.fid(1,1) = sp.fid(1,1) .* 2;

scaleFMM = 1.8482; % Gosia calibration scale factor
refppm = 4.66;

% Shift both
sp.ppmref = refppm;
spW.ppmref = refppm;

% Uncomment these to look at the data
%view_spectrum_gui(sp);
%view_spectrum_gui(spW);


%% Read basis set
% these have TSP at  0, so shift is 4.66
[spBasis, labels] = readMRSChallengeBasisSet();
spBasis.ppmref = refppm;

%% Need a simulated water for fitting water
llw = 3; % Looks like the basis sets are 3Hz
glw = 0;
spref = spW;
spref.ppmref = 0;
[~, water] = modelMultipleVoigts([1, 0, 0, llw, glw], [], spref);
water.ppmref = refppm;

% Step 1: Fit the water

% Scale it up so the fit parameters are ok
scaleW = 10000;

% Setup a fit
theta0 = thetaInitStruct(1,1);

theta0.global(2).value = 0;
theta0.region(1).element(1).amp = 1;
theta0.region(1).element(1).phs = 0;
theta0.region(1).element(1).llw = 10;
theta0.region(1).element(1).glw = 10;
theta0.region(1).element(1).delf = 0;

theta0Arr = thetaStructToArray(theta0);
thetaPrint(theta0);

residual = lcSpecResidual_v3(theta0Arr, spW, water.*scaleW);

%%
if bPrePhaseWater
   [spW phsrad0] = spW.aph0;
   fprintf('Pre-phasing water found %.1f degrees\n', phsrad0*180/pi);
end

[thetaArrW, ~, spSynW, ncrbW] = fitSpecLC_v3water( spW, theta0Arr, water .* scaleW);
spSynW = spSynW .* scaleW;
thetaW = thetaArrayToStruct(thetaArrW, theta0);
thetaPrint(thetaW);

Awater = thetaW.region(1).element(1).amp * scaleW;
phsWater = thetaW.region(1).element(1).phs;
llwWater = thetaW.region(1).element(1).llw;
glwWater = thetaW.region(1).element(1).glw;
fprintf('\n*** Results: amp=%.3f, phs=%.1f deg, llw=%.1f Hz, glw=%.1f Hz ***\n', ...
    Awater, phsWater, llwWater, glwWater);


%% Plot the water
figure(9)
spdiffW = spW.*scaleW - spSynW;
resOff = -0;
plot(spW.freqPPM, real(spW.spec), '-k', ...
    spSynW.freqPPM, real(spSynW.spec), '-b', ...
    spdiffW.freqPPM, real(spdiffW.spec)+resOff, '-r');
legend('data', 'fit', 'residual')
set(gca, 'xdir', 'reverse');
set(gca, 'xlim', [4.3 5]);

set(gca, 'yticklabel', [])
xlabel('ppm')
title(sprintf('dataset #%d Water (%.1f/%.1f Hz)', dsNum, llw, glw))

saveas(gcf, fnameWaterFit, 'jpg');

%% Now fit metabolite spectrum

% Optionally add Lipid ot the basis set
if bIncludeLipid
    labels{size(labels,2)+1} = 'Lip';
    
    % This basis set I've derived from the measured data.
    cf = 4.66;
    llw = 10;
    glw = 10;
    [~, lipid] = modelMultipleVoigts(...
        [0.642, 1.30-cf, 0, llw, glw, ...
        0.088, 0.90-cf, 0, llw, glw, ...
        0.062, 2.02-cf, 0, llw, glw, ...
        0.058, 1.60-cf, 0, llw, glw, ...
        0.058, 2.24-cf, 0, llw, glw, ...
        0.037, 5.29-cf, 0, llw, glw, ...
        0.010, 5.19-cf, 0, llw, glw, ...
        0.006, 2.75-cf, 0, llw, glw, ...
        0.039, 4.20-cf, 0, llw, glw, ...
        0.005, 3.20-cf, 0, llw, glw*2], ...
        [], spref);
    %lipid = lipid.shiftFreq(refppm);
    spBasisMod = spBasis.append(lipid);
else
    spBasisMod = spBasis;
end



%%
% Setup a fit
theta0 = thetaInitStruct(1,spBasisMod.M);
for idx=1:spBasisMod.M; theta0.region(1).element(idx).label = labels{idx};end;

theta0.global(1).value = phsWater;
theta0.global(2).value = llwWater;
theta0.global(3).value = glwWater;


theta0.region(1).element(1).amp = 1;
theta0.region(1).element(2).amp = 1;
theta0.region(1).element(3).amp = 1;
theta0.region(1).element(4).amp = 3;
theta0.region(1).element(5).amp = 1;
theta0.region(1).element(6).amp = 1;
theta0.region(1).element(7).amp = 1;
theta0.region(1).element(8).amp = 1;
theta0.region(1).element(9).amp = 1;
theta0.region(1).element(10).amp = 1;
theta0.region(1).element(11).amp = 1;
theta0.region(1).element(12).amp = 1;
theta0.region(1).element(13).amp = 0;
theta0.region(1).element(14).amp = 1;
theta0.region(1).element(15).amp = 8;
theta0.region(1).element(16).amp = 15;
theta0.region(1).element(17).amp = 1;
theta0.region(1).element(18).amp = 8;
theta0.region(1).element(19).amp = 1;
theta0.region(1).element(20).amp = 1;
theta0.region(1).element(21).amp = 1;
if bIncludeLipid
    % Then I ramp up both lipid and lacatate
    theta0.region(1).element(13).amp = 100;
    theta0.region(1).element(22).amp = 10000;
end

thetaPrint(theta0);
theta0Arr = thetaStructToArray(theta0);


%% HSVD editing of residual water
if bHSVDWaterRemoval
    %[residual, peaks, spfit] = hsvdRemoveInside(sp, [4.2 5.2], 1000);
    [residual, peaks, spfit] = hsvdRemoveInside(sp, [3.9 4.5], 1000);
    sp = residual;
    
    % Lets compare
    figure(21)
    plot(sp.freqPPM, real(sp.spec), '-k', ...
        residual.freqPPM, real(residual.spec), '-b');
    legend('Original', 'edited');
    set(gca, 'xlim', [0 6]);
    set(gca, 'xdir', 'reverse');
    title(sprintf('HSVD residual water removal (#%d)', dsNum));
end

%% Prephase with water phase
if bPrePhaseWater
   sp = sp.phaseShift(-phsrad0); 
end

%%
[sthetaArr, ~, spSyn, ncrbM] = fitSpecLC_v3metabolite( sp, theta0Arr, spBasisMod);

phsGlb = sthetaArr(1);
llwGlb = sthetaArr(2);
glwGlb = sthetaArr(3);
delfGlb  = sthetaArr(4);
amps = sthetaArr(5:end);
fprintf('\nGlobals: phs=%.1f, llw=%.1f, glw=%.1f, defl=%.1f\n', phsGlb, llwGlb, glwGlb, delfGlb);

%% Quantify
WConcCorr = 29697;
MetCorr = 1.206;
totalScaleF = scaleFMM * WConcCorr * MetCorr / Awater;
fprintf('\n*** Results ***\n')
fp = fopen(fnameResults, 'w');
for idx=1:spBasisMod.M
    %conc = thetaMet.region(1).element(idx).amp * totalScaleF;
    conc = amps(idx) * totalScaleF;
    %fprintf('%d [%s]: \t%.3f mM\n', idx, labels{idx}, conc);
    fprintf('%s\t%.3f\n', labels{idx}, conc);
    fprintf(fp, '%s\t%.3f\n', labels{idx}, conc); 
end
fprintf('\nResults also written to %s\n', fnameResults);
fclose(fp);

%% Display
figure(10)
spdiff = sp - spSyn;
resOff = -0;
plot(sp.freqPPM, real(sp.spec), '-k', ...
    spSyn.freqPPM, real(spSyn.spec), '-b', ...
    spdiff.freqPPM, real(spdiff.spec)+resOff, '-r');
legend('data', 'fit', 'residual')
set(gca, 'xdir', 'reverse');
set(gca, 'xlim', [0 4.3]);
%set(gca, 'xlim', [-2 7]);
set(gca, 'yticklabel', [])
xlabel('ppm')
title(sprintf('dataset #%d (%.1f/%.1f Hz)', dsNum, llwGlb, glwGlb))

saveas(gcf, fnameMetaboliteFit, 'jpg');














