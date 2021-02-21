% Function that performs NLLS fitting of the an MBSpectrum using a basis set 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version does not use frequency selectivity

% v3 is v2 with Voigt lineshapes
%   there is a water and a metabolite version. Different fitting ranges and
%   parameter limits
% v2 is tailored for the MRS fitting contest


function [stheta, resnorm, spSyn, ncrb] = fitSpecLC_v3metabolite( spdata, theta0, spbasis)

assert(spdata.N == spbasis.N);
assert(spdata.frq == spbasis.frq);
%assert(spdata.ppmref == spbasis.ppmref);
assert(spdata.M == 1);

J = spbasis.M;
nGlobals = 4;
nParamsPerElement = 5;
assert(J == (size(theta0,2)-nGlobals)/nParamsPerElement);

% As a hack I'm going to do a "small" theta which is ONLY amplitudes per
% element
stheta0 = zeros(J,1);
lb = stheta0 .* 0;
ub = stheta0 .* 0;
for jdx = 1:J
    
    lb(jdx) = 0;
    ub(jdx) = Inf;

    % joff is the offset into the larger theta0 array
    joff = nGlobals + nParamsPerElement*(jdx-1);
    stheta0(jdx) = theta0(joff+1);
end

% Prepand Globals
phsGlb = theta0(1);
llwGlb = theta0(2);
glwGlb = theta0(3);
delfGlb = 0;
stheta0 = [phsGlb; llwGlb; glwGlb; delfGlb; stheta0];
lb = [-2; 0; 0; -50; lb];
ub = [2; 10; 10; 50; ub];


% Make an anonymous function to pass additional parameters
rangePPM = [0.5 4.3];
anonFun = @(x)lcSpecResidual_v3metabolite(x, spdata, spbasis, rangePPM); 


% Specify fitting options
options = optimoptions(@lsqnonlin,  'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter-detailed', 'TolFun', 1E-8, 'TolX', 1E-8, ...
    'MaxFunEvals', 20000);

% Begin Fitting
start = tic;
[stheta, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(anonFun, stheta0, lb, ub, options);
dur = toc(start);
fprintf('Completed in %.1f s\n', dur);


disp(output)

% % Variance and Cramer Rao
% % I'm only going to do this on amplitudes
% [~, ~, noiserms, ~] = getSnrFd(spdata);
% noiserms_real = noiserms / sqrt(2);
% varb = full(inv((1/noiserms_real^2)* jacobian' * jacobian));
% crb = sqrt(diag(varb));
% for idx=1:spbasis.M
%     index = 3 + 4*(idx-1) + 1;
%     ncrb(idx) = crb(index,1) / abs(theta(index)); 
%     %ncrb(idx) = crb(index,1) / abs(theta(index)); 
% end
ncrb = [];

% Synthesize the fit
[~, spSyn] = lcSpecResidual_v3metabolite(stheta, spdata, spbasis, rangePPM);

fprintf('Resnorm = %.1f, Resnorm per pt = %.1f\n', resnorm, resnorm/numel(residual));

% Calculate the data norm
tmp = [real(spdata.spec); imag(spdata.spec)];
datanorm = sum(tmp.^2);
relresnorm = resnorm/datanorm;
fprintf('Datanorm = %.1f, RelRes = %.2f%%\n', datanorm, relresnorm*100);





