% Function that performs NLLS fitting of the an MBSpectrum using a basis set 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First version does not use frequency selectivity, has flexibility geared
% toward my SLIM solutions

% v3 is v2 with Voigt lineshapes
%   there is a water and a metabolite version. Different fitting ranges and
%   parameter limits
% v2 is tailored for the MRS fitting contest


function [theta, resnorm, spSyn, ncrb] = fitSpecLC_v3water( spdata, theta0, spbasis)

assert(spdata.N == spbasis.N);
assert(spdata.frq == spbasis.frq);
%assert(spdata.ppmref == spbasis.ppmref);
assert(spdata.M == 1);

J = spbasis.M;
nGlobals = 4;
nParamsPerElement = 5;
assert(J == (size(theta0,2)-nGlobals)/nParamsPerElement);

% Bounds. Must match theta0 structure
lb = theta0 .* 0;
ub = theta0 .* 0;

% YOu can effectively disable a parameter by giving it a tiny range

% ranges
phsMax = .0001; % Degrees
llwMin = 0; 
llwMax = .0001;
glwMin = 0;
glwMax = .0001;
delfMax = .0001; % Hz


% Elemental terms
phsEMax = 179;
llwEMin = 0; 
llwEMax = 20;
glwEMin = 0;
glwEMax = 20;
delfEMax = 20;

% Globals
lb(1) = -phsMax;
ub(1) = phsMax;

lb(2) = llwMin;
ub(2) = llwMax;

lb(3) = glwMin;
ub(3) = glwMax;

lb(4) = -delfMax;
ub(4) = delfMax;

% Per basis component j
for jdx=1:J
    joff = nGlobals + nParamsPerElement*(jdx-1);
    
    % amp
    lb(joff + 1) = 0;
    ub(joff + 1) = Inf;
    
    % phs
    lb(joff + 2) = -phsEMax;
    ub(joff + 2) = phsEMax;
        
    % llw
    lb(joff + 3) = llwEMin;
    ub(joff + 3) = llwEMax;

    % glw
    lb(joff + 4) = glwEMin;
    ub(joff + 4) = glwEMax;
    
    % delf
    lb(joff + 5) = -delfEMax;
    ub(joff + 5) = delfEMax;
    
end


% Make an anonymous function to pass additional parameters
rangePPM = [4.4 5.0];
anonFun = @(x)lcSpecResidual_v3(x, spdata, spbasis, rangePPM); 


% Specify fitting options
options = optimoptions(@lsqnonlin,  'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter-detailed', 'TolFun', 1E-6, 'TolX', 1E-6, ...
    'MaxFunEvals', 200);

% Begin Fitting
start = tic;
[theta, resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(anonFun, theta0, lb, ub, options);
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
[~, spSyn] = lcSpecResidual_v3(theta, spdata, spbasis);

fprintf('Resnorm = %.1f, Resnorm per pt = %.1f\n', resnorm, resnorm/numel(residual));

% Calculate the data norm
tmp = [real(spdata.spec); imag(spdata.spec)];
datanorm = sum(tmp.^2);
relresnorm = resnorm/datanorm;
fprintf('Datanorm = %.1f, RelRes = %.2f%%\n', datanorm, relresnorm*100);





