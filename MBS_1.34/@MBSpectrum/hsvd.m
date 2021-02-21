% hsvd - Performs a Henkel Singular Value Decomposition of a spectrum
%
% This is a very basic HSVD decomposition as described by Laudadio et
% al, JMR 2002. Note that is not the original HSVD reference but a clearly
% described implementation.
%
%   [residual, peaks, spfit] = hsvd(sp, nMaxValues)
%       sp - input MBSpectrum
%       nMaxVales - maximum values to decompose into. Optional
%       residual - an MBSpectrum of the residual
%       peaks - an array of lorentzian peak values. See simulateLorentzian
%           for structure
%       spfit - an NxM MBSpectrum of all M singular values
%
%   Note: This only works with a single (M=1) spectrum. 
%   Note: The performance slows down by N^2. Truncate spe

% Created: 8/6/2010 Patrick Bolan and Timo Liimatainen
% MBS - Minnesota Breast Spectroscopy package
function [residual, peaks, spfit] = hsvd(sp, nMaxValues)

if nargin<2
    nMaxValues = Inf;
end

% Needs to operate on 1 at a time
if (sp.M>1)
    error('hsvd only operates on 1 at a time');
end

% Extract the FID once
y = sp.fid.';

% Make a Hankel matrix, taking the upper left subspace.
% Eq [2.2]
L = sp.N/2;
M = sp.N-L-1;
H = zeros(L, M);
for idx = 1:L
    H(idx, :) = y(idx + (1:M));
end


% SVD, Eq [2.3]
fprintf('Starting svd...');
[~, ~, V]=svd(H);
fprintf(' done.\n');

% Need to pick rank K. Should be 1 more than the number of peaks I want.
% The bigger this is the better, and M is max
K = M;

Vk = V(1:K, :);
Vb = Vk(2:K, :);
Vt = Vk(1:K-1, :);

% Calc least-squares soln
E = Vt \ Vb;

% Step 4
zk = eig(E);
zk = zk(1:K-1);

% Step 5
% Must be careful with this atan2 function
fk = (atan2(real(zk), imag(zk)) - pi/2) / (2 * pi * sp.dt);
dk = -log(abs(zk)) / sp.dt;

% eliminate those with negative d values 
goodIndices = dk>0;
fk = fk(goodIndices);
dk = dk(goodIndices);

%fk = -1 .* fk;

% Now to solve for ck = ak * exp(i*phik), the complex amplitude of each
% singular value, we make the eigenvector M and solve.  Timo called this M
% matrix 'H', I don't know why.
M = zeros(length(dk), length(sp.time));
for idx=1:length(dk)
   %M(idx,:) = exp((-dk(idx) + 1i*2*pi*fk(idx)) .* sp.time); 
   M(idx,:) = exp((-dk(idx) + 1i*2*pi*fk(idx)) .* sp.time); 
end

ck = M.' \ y.';
% ck = pinv(M.')*y.'; % Same thing, maybe more robust
ak = abs(ck);
phik = imag(log(ck./ak));


% Sort these by amplitude
[ak, sortedIdx] = sort(ak, 'descend');
fk = fk(sortedIdx);
phik = phik(sortedIdx);
dk = dk(sortedIdx);

% Select a subset of values for calculating residuals
nVals = min(size(ak,1), nMaxValues);
ak = ak(1:nVals);
fk = fk(1:nVals);
phik = phik(1:nVals);
dk = dk(1:nVals);


% Prepare output structure
for idx=1:nVals
   peaks(idx).A = ak(idx);
   %peaks(idx).lw = 1./(d(idx)*pi);
   %peaks(idx).lambda = pi./(dk(idx));
   peaks(idx).lambda = dk(idx); % Lambda = FWHM*pi
   peaks(idx).freq = -fk(idx); % Note the negative sign here!
   peaks(idx).phs = phik(idx);
end

% Report the values found
% fprintf('Found %d acceptable singular values:\n', nVals);
% for idx=1:nVals
%     fprintf('   %d: A = %8.2f \tfreq = %.4f (Hz) \tlambda = %.4f (1/s)\tphi = %.4f (deg)\n', ...
%        idx, peaks(idx).A, peaks(idx).freq, peaks(idx).lambda, peaks(idx).phs*180/pi);
% end

% Create MBSpectra from these peaks
spfit = simulateLorentzian(peaks, sp.N, sp.sw, sp.frq);
residual = sp - spfit.sum;


% % For reference, here's some diagnostic plots
% figure(10)
% plotSpecArray(spfit);
% 
% figure(11)
% freq = sp.freqHz;
% plot(freq, real(sp.spec), '-b', ...
%     freq, real(spfit.sum.spec), '-k', ...
%     freq, real(residual.spec), '-r');
% set(gca, 'xdir', 'reverse');

return;


% % Simulate the solution
% freq = sp.freqHz;
% figure(9)
% clf
% title('Components')
% hold on
% fitfid = zeros(nVals, sp.N);
% for idx=1:nVals
%    fitfid(idx,:) = ak(idx) * exp(1i * phik(idx)) .* exp(sp.time.*(-dk(idx) + 1i*2*pi * fk(idx))); 
%    plot(freq, real(fftshift(fft(fitfid(idx,:)))));
% end
% hold off
% 
% fullfit = sum(fitfid,1);
% residual = y - fullfit;
% figure(10)
% plot(freq, real(fftshift(fft(y))) ./ sqrt(sp.N), '-b', ...
%     freq, real(fftshift(fft(fullfit))) ./ sqrt(sp.N), '-k', ...
%     freq, real(fftshift(fft(residual))) ./ sqrt(sp.N), '-r');
% set(gca, 'xdir', 'reverse');
% legend('data', 'fit', 'residual');






