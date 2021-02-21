% ******************************************************************************
%             MBS - Minnesota Breast Spectroscopy analysis package
%               Developed by Patrick Bolan and Michael Garwood
% ******************************************************************************
% FUNCTION: mbsSpectrum - leftTrimPhilipsMaxEcho
% AUTHOR: pjb
% CREATED: 2/17/2009
% DESCRIPTION: Philips has a bug in their current press sequence in that
% the receiver is gated on at the top of the echo. This is a problem
% because it takes 2-5 pts for the receiver to stabilize. One solution is
% to start acquiring after the last gradient, and then throw away points
% prior to the theoretical top of the echo. Philips has a MAXECHO setting
% that lets you do just that.
%
% In the SDAT file, t0_mu1_direction is a percentage of points that are
% left of the nominal echo. So you should just be able to figure out how
% many points to throw away and you're done! But, the # of points is not
% necesicarily integer, and the linear phase gets messed
%
% ARGUMENTS: a philips mbsSpectrum
% RETURNS: mbsSpectrum
% MODIFICATIONS:
% ******************************************************************************
function out = leftTrimPhilipsMaxEcho(sp)

% The "max echo" setting was used. Throw out some samples to the left,
% and append zeros to the right.

% I'll do this using the shift theorem: apply a linear phase to the
% spectrum and that'll shift the fid. This allows sub-pixel shifts,
% effectively with a sinc interpolation.

% Note this also changes the spectrum size to something other than 2^n,
% which then changes your FFTs into DFTs (slow!).

% t0_mu1_direction is a percentage of points that are left of the
% nominal echo.
%num_pts_to_clip = floor(sp.header.t0_mu1_direction * sp.header.samples/100);
num_pts_to_clip = sp.header.t0_mu1_direction * sp.header.samples/100;
fprintf('clipping %.3f points from left of FID\n', num_pts_to_clip);

% You need to do a DC correction first or this can make gibbs ringing
sp = correctDC(sp);

% This shifts left, sub-point, and makes the whole thing smaller
% I call this just to figure out how many point will be clipped.
tmp = pretruncate(sp, num_pts_to_clip);

% Trying to avoid the wrap artifact. 
origpts = sp.N;
finalpts = tmp.N;


% HACK! This is a fudge factor based on phantom scans. Improves the phase
% slightly 
num_pts_to_clip = num_pts_to_clip * 1.002;


sp_zf2 = zeroFill(sp, 2*origpts);
shifted_zf2 = pretruncate(sp_zf2, num_pts_to_clip);
out = zeroFill(shifted_zf2, finalpts);

% Or, should it be the same as the original? Still figuring out what is
% best
out = zeroFill(out, origpts);

% The shifting is implemented in shiftFid()
% sp = shiftFid(sp, -num_pts_to_clip);
% 
% % Now truncate. The new size should be >= num_pts_to_clip and even
% newsize = sp.pts - ceil(num_pts_to_clip/2)*2; 
% 
% disp(sprintf('t0_mu1_direction = %.4f, shift %.4f out of %d points, new size %d',...
%     sp.header.t0_mu1_direction, num_pts_to_clip, sp.header.samples, newsize));
% 
% sp = truncate(sp, newsize);
return;


% % PJB 20090213
% % When the t0_mu1_direction was getting large, it didn't look quite
% % right. It was throwing away the highest point, and some linear phase
% % was left in the spectrum.
% 
% warning('MBS:PhilipsMaxEcho', ...
%     'MAXECHO setting used, discarding %.0d points on left', num_pts_to_clip);
% 
% % This is the version with zero-filling on the right. Produces artifacts
% %tmp = sp.fid .* 0;
% %tmp(1:sp.pts-num_pts_to_clip, :) = sp.fid(num_pts_to_clip+1:sp.pts, :);
% %tmp(sp.pts-num_pts_to_clip+1:sp.pts, :) = 0;
% 
% % Here is the non-zero-filling version. The spectrum shrinks
% tmp = sp.fid(num_pts_to_clip+1:sp.pts, :);
% sp.fid = tmp;
% pts = size(tmp,1);
% 
% % Adjust header info
% sp.at = sp.at * pts/sp.pts;
% sp.pts = pts;
% 
% % Regenerate spec and axes
% sp = spec_fft(sp);
% sp = calcAxes(sp);
